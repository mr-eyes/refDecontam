#include "string"
#include "kDataFrame.hpp"
#include <stdexcept>
#include "tuple"
#include <sys/stat.h>
#include "colored_kDataFrame.hpp"
#include <zlib.h>
#include <cstdio>
#include <kseq/kseq.h>
#include <iostream>
#include <vector>
#include <cassert>
#include "kmerDecoder.hpp"

#define KSIZE 21
#define CHUNK_SIZE 5000

using namespace std;
using namespace phmap;


string create_dir(string output_file, int serial) {
    int dir_err;
    string new_name = "";

    if (!serial) {
        dir_err = mkdir(output_file.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        new_name = output_file;
    }
    else {
        new_name = output_file + "_v." + to_string(serial);
        dir_err = mkdir(new_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    if (-1 == dir_err) return create_dir(output_file, ++serial);

    return new_name;
}

class fileHandler {

public:
    ofstream fileStream;

    fileHandler(string& filename) {
        this->fileStream.open(filename);
    }

    void write(string& line) {
        this->fileStream << line;
    }

    void close() {
        fileStream.close();
    }

};

inline bool file_exists(const std::string& name) {
    struct stat buffer
    {
    };
    return (stat(name.c_str(), &buffer) == 0);
}

inline string time_diff(std::chrono::high_resolution_clock::time_point& t1) {
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    long hr = milli / 3600000;
    milli = milli - 3600000 * hr;
    long min = milli / 60000;
    milli = milli - 60000 * min;
    long sec = milli / 1000;
    milli = milli - 1000 * sec;
    string timeDiff;
    timeDiff.append(to_string(min));
    timeDiff.append(":");
    timeDiff.append(to_string(sec));
    timeDiff.append(":");
    timeDiff.append(to_string(milli));

    return timeDiff;
}

flat_hash_map<uint64_t, std::vector<uint32_t>> load_colors(string index_prefix) {
    flat_hash_map<uint64_t, std::vector<uint32_t> > colors;
    string inputFilename = index_prefix + "colors.intvectors";
    ifstream input(inputFilename);
    uint32_t size;
    input >> size;
    colors = flat_hash_map<uint64_t, std::vector<uint32_t> >(size);
    for (int i = 0; i < size; i++) {
        uint64_t color, colorSize;
        input >> color >> colorSize;
        uint32_t sampleID;
        colors[color] = std::vector<uint32_t>(colorSize);
        for (int j = 0; j < colorSize; j++) {
            input >> sampleID;
            colors[color][j] = sampleID;
        }
    }
    return colors;
}

tuple<string, vector<int>> score(vector<uint32_t>& genomes) {

    vector<int> sources;

    if (genomes.empty())
        return make_tuple("unmapped", sources);

    flat_hash_map<int, int> scores;
    flat_hash_map<int, int> reverse_scores;
    flat_hash_map<int, int> countFreq;
    vector<int> all_scores;

    for (const auto& genome : genomes) {
        scores[genome]++;
    }

    for (const auto& score : scores) {
        countFreq[score.second]++;
        all_scores.emplace_back(score.second);
        reverse_scores[score.second] = score.first;
    }

    auto max = std::max_element(all_scores.begin(), all_scores.end());

    if (countFreq[*max] == 1) {
        sources.emplace_back(reverse_scores[*max]);
        return make_tuple("unique", sources);
    }

    for (const auto& score : scores) {
        if (score.second == *max) {
            sources.emplace_back(score.first);
        }
    }

    return make_tuple("ambig", sources);

}



int main(int argc, char** argv) {

    KS_FULL_COMMENT = true;

    if (argc != 3) {
        cerr << "run ./decontaminate <ref_fasta> <index_prefix>" << endl;
        exit(1);
    }

    string reads_file = argv[1];
    string index_prefix = argv[2];

    if (!file_exists(reads_file)) {
        throw std::runtime_error("Could not open the unitigs fasta file");
    }

    if (!file_exists(index_prefix + ".extra")) {
        throw std::runtime_error("Could not open kProcessor index file");
    }

    // kProcessor Index Loading
    std::cerr << "Loading kProcessor index ..." << std::endl;
    colored_kDataFrame* ckf = colored_kDataFrame::load(index_prefix);
    kDataFrame* kf = ckf->getkDataFrame();

    set<int> vec_singleColors;
    flat_hash_map<uint64_t, vector<uint32_t>> color_to_vecGroups;
    flat_hash_map<uint64_t, string> color_to_groupString;
    cerr << "Loading colors .." << endl;
    auto colorsIntVector = load_colors(index_prefix);
    for (auto const& color : colorsIntVector) {
        uint64_t color_id = color.first;
        auto all_group_ids = color.second;

        for (auto _grp_id : all_group_ids) {
            color_to_vecGroups[color_id].emplace_back(_grp_id);
            vec_singleColors.emplace(_grp_id);
        }
    }

    map<string, fileHandler*> fasta_writer;

    cerr << "Creating fasta file handlers..." << endl;
    for (auto& item : vec_singleColors) {
        string file_name = "genome_" + ckf->namesMap[item] + "_partition.fa";
        fasta_writer[to_string(item)] = new fileHandler(file_name);
    }

    string unmapped_file_name = "unmapped_partition.fa";
    fasta_writer["unmapped"] = new fileHandler(unmapped_file_name);

    int chunkSize = CHUNK_SIZE;
    int kSize = KSIZE;


    gzFile fp;
    kseq_t* kseqObj;
    fp = gzopen(reads_file.c_str(), "r");
    kseqObj = kseq_init(fp);

    cout << "Processing started ..." << endl;

    auto * KD = kmerDecoder::getInstance(KMERS, mumur_hasher, {{"kSize", kSize}});
    

    int total = 0;
    int chunks = 0;

    uint64_t _unique = 0;
    uint64_t _ambig = 0;
    uint64_t _unmatched = 0;

    while (kseq_read(kseqObj) >= 0) {
        // if (string(kseqObj->seq.s).size() < kSize) continue;

        std::string seq = kseqObj->seq.s;
        std::string id;
        id.append(kseqObj->name.s);
        if (kseqObj->comment.l) id.append(kseqObj->comment.s);

        std::string record = ">";
        record.append(id);
        record.append("\n");
        record.append(seq);
        record.append("\n");

        vector<uint32_t> kmers_matches;

        for (unsigned long i = 0; i < seq.size() - kSize + 1; i++) {
            uint64_t kmer = KD->hash_kmer(seq.substr(i, kSize));
            uint64_t color = kf->getCount(kmer);
            for (const auto& genomeID : color_to_vecGroups[color]) {
                kmers_matches.emplace_back(genomeID);
            }
        }

        auto category = score(kmers_matches);


        if (get<0>(category) == "unique") {
            _unique++;
            fasta_writer[to_string(get<1>(category)[0])]->write(record);
        }        
else if (get<0>(category) == "unmapped") {
            _unmatched++;
            fasta_writer["unmapped"]->write(record);
        }
        else if (get<0>(category) == "ambig") {
            _ambig++;
            for (auto const& genomeID : get<1>(category)) {
                fasta_writer[to_string(genomeID)]->write(record);
            }
        }


        total++;
        if (total == 5000) {
            cout << "processed: " << 5000 * ++chunks << "contigs" << endl;
            total = 0;
        }

    }

    cout << "_unique: " << _unique << endl;
    cout << "_ambig: " << _ambig << endl;
    cout << "_unmatched: " << _unmatched << endl;

    for (auto f : fasta_writer)
        f.second->close();


    return 0;
}