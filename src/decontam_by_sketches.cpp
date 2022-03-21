#include "kDataFrame.hpp"
#include <stdexcept>
#include "tuple"
#include <sys/stat.h>
#include "colored_kDataFrame.hpp"
#include <kseq/kseq.h>
#include <zlib.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <glob.h> // glob(), globfree()
#include <cassert>
#include "kmerDecoder.hpp"
#include <string>     // std::string, std::to_string

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
        new_name = output_file + "_v." + std::to_string(serial);
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

// thanks to https://stackoverflow.com/a/8615450/3371177
std::vector<std::string> glob(const std::string& pattern) {
    using namespace std;

    // glob struct resides on the stack
    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    // do the glob operation
    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if (return_value != 0) {
        globfree(&glob_result);
        stringstream ss;
        ss << "glob() failed with return_value " << return_value << endl;
        throw std::runtime_error(ss.str());
    }

    // collect all the filenames into a std::list<std::string>
    vector<string> filenames;
    for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
        filenames.push_back(string(glob_result.gl_pathv[i]));
    }

    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}

tuple<string, vector<int>> score(vector<uint32_t>& genomes) {

    vector<int> sources;

    if (genomes.empty())
        return make_tuple("unmapped", sources);

    flat_hash_map<int, int> genome_to_count;
    flat_hash_map<int, int> count_to_genome;
    flat_hash_map<int, int> countFreq;
    vector<int> all_counts;

    // Just count frequencies
    for (const auto& genome : genomes) {
        genome_to_count[genome]++;
    }

    // Count max freq
    for (const auto& item : genome_to_count) {
        countFreq[item.second]++;
        all_counts.emplace_back(item.second);
        count_to_genome[item.second] = item.first;
    }

    auto max_count = std::max_element(all_counts.begin(), all_counts.end());

    // All kmers are coming from a single genome
    if (countFreq[*max_count] == 1) {
        sources.emplace_back(count_to_genome[*max_count]);
        return make_tuple("unique", sources);
    }

    // Kmers are coming from different genomes
    for (const auto& score : genome_to_count)
        if (score.second == *max_count)
            sources.emplace_back(score.first);

    return make_tuple("ambig", sources);
}


tuple<vector<int>, int> score_by_max(vector<uint32_t>& genomes) {

    vector<int> sources;

    if (genomes.empty()) return make_tuple(sources, 0);

    // Count frequencies
    flat_hash_map<uint32_t, uint32_t> genome_to_freq;
    for (auto const& item : genomes) genome_to_freq[item]++;

    // Only max genomes
    vector<int> freqs;

    for (auto& it : genome_to_freq) freqs.emplace_back(it.second);

    // Max freq
    auto max = std::max_element(freqs.begin(), freqs.end());

    // Get all genomes with freq = max
    for (auto& item : genome_to_freq)
        if (item.second == *max)
            sources.emplace_back(item.first);

    return make_tuple(sources, *max);
}


int main(int argc, char** argv) {

    KS_FULL_COMMENT = true;

    if (argc != 3) {
        cerr << "run ./decontam_by_sketches <ref_fasta> <sketches_dir>" << endl;
        exit(1);
    }

    string reads_file = argv[1];
    string kfs_dir = argv[2];

    flat_hash_map<string, flat_hash_map<string, int>> stats;

    kDataFrame* frame;
    std::string dir_prefix = kfs_dir.substr(kfs_dir.find_last_of("/\\") + 1);

    flat_hash_map<string, string> namesMap;
    string names_fileName = kfs_dir;

    string line;

    int total_kfs_number = 0;

    // get kSize and type
    for (const auto& dirEntry : glob(kfs_dir + "/*")) {
        string file_name = (string)dirEntry;
        size_t lastindex = file_name.find_last_of(".");
        string kf_prefix = file_name.substr(0, lastindex);
        std::string::size_type idx;
        idx = file_name.rfind('.');
        std::string extension = "";
        if (idx != std::string::npos) extension = file_name.substr(idx + 1);
        int detected_kSize;
        hashingModes _hm;
        if (extension == "mqf" || extension == "phmap") {
            auto* _kf = kDataFrame::load(kf_prefix);
            _hm = _kf->getkmerDecoder()->hash_mode;
            detected_kSize = _kf->getkSize();
            cout << "Detected kSize: " << detected_kSize << endl;
        }
        else {
            continue;
        }
        // if(extension == "mqf") {frame = new kDataFrameMQF(detected_kSize, 30, mumur_hasher); break;} // temp. switch off
        if (extension == "mqf") { frame = new kDataFramePHMAP(detected_kSize, _hm); break; }
        else if (extension == "phmap") { frame = new kDataFramePHMAP(detected_kSize, _hm); break; }
        else { continue; }
    }


    if (!file_exists(reads_file)) {
        throw std::runtime_error("Could not open the reads fasta file");
    }


    flat_hash_map<string, kDataFrame*> all_kfs;
    flat_hash_map<string, int> genome_to_id;
    flat_hash_map<int, string> id_to_genome;

    // Loading kfs
    int counter = 0;
    for (const auto& dirEntry : glob(kfs_dir + "/*")) {
        string file_name = (string)dirEntry;
        size_t lastindex = file_name.find_last_of(".");
        string kf_prefix = file_name.substr(0, lastindex);
        std::string kf_basename = kf_prefix.substr(kf_prefix.find_last_of("/\\") + 1);
        std::string::size_type idx;
        idx = file_name.rfind('.');
        std::string extension = "";
        if (idx != std::string::npos) extension = file_name.substr(idx + 1);
        if (extension != "mqf" and extension != "phmap") continue;
        cout << "loading " << kf_prefix << "..." << endl;
        all_kfs[kf_basename] = kDataFrame::load(kf_prefix);
        genome_to_id[kf_basename] = counter++;
        id_to_genome[genome_to_id[kf_basename]] = kf_basename;
        stats[kf_basename] = { {"unique", 0}, {"ambig", 0} };
    }

    map<string, fileHandler*> fasta_writer;

    cerr << "Creating fasta file handlers..." << endl;
    for (auto& item : all_kfs) {
        string file_name = "genome_" + item.first + "_partition.fa";
        fasta_writer[item.first] = new fileHandler(file_name);
    }

    string unmapped_file_name = "unmapped_partition.fa";
    fasta_writer["unmapped"] = new fileHandler(unmapped_file_name);
    string ambig_file_name = "ambiguous_partition.fa";
    fasta_writer["ambig"] = new fileHandler(ambig_file_name);

    int chunkSize = CHUNK_SIZE;
    int kSize = KSIZE;


    gzFile fp;
    kseq_t* kseqObj;
    fp = gzopen(reads_file.c_str(), "r");
    kseqObj = kseq_init(fp);

    cout << "Processing started ..." << endl;

    auto* KD = kmerDecoder::getInstance(KMERS, mumur_hasher, { {"kSize", kSize} });


    int total = 0;
    int chunks = 0;

    uint64_t _unmatched = 0;

    while (kseq_read(kseqObj) >= 0) {
        // if (string(kseqObj->seq.s).size() < kSize) continue;

        std::string seq = kseqObj->seq.s;
        std::string id;
        id.append(kseqObj->name.s);
        if (kseqObj->comment.l) id.append(kseqObj->comment.s);

        vector<uint32_t> kmers_matches;

        phmap::flat_hash_set<uint64_t> hashes_set;
        

        // unique hashes once
        for (unsigned long i = 0; i < seq.size() - kSize + 1; i++)
            hashes_set.insert(KD->hash_kmer(seq.substr(i, kSize)));
        
        // Iterate and search over all kfs
        for (auto const& ref : all_kfs)
            for (const auto& kmer_hash : hashes_set)
                if (ref.second->getCount(kmer_hash))
                    kmers_matches.emplace_back(genome_to_id[ref.first]);

        tuple<vector<int>, int> scored_genomes = score_by_max(kmers_matches);
        if (get<0>(scored_genomes).size() == 1) {
            double percentage = get<1>(scored_genomes) / (double)hashes_set.size();
            string header_tail = "|";
            header_tail.append("kmers:" + to_string(get<1>(scored_genomes)) + ";");
            header_tail.append("percentage:" + to_string(percentage));
            std::string record = ">";
            // TODO: Refactor later.
            record.append(id);
            record.append(header_tail);
            record.append("\n");
            record.append(seq);
            record.append("\n");
            stats[id_to_genome[get<0>(scored_genomes)[0]]]["unique"]++;
            fasta_writer[id_to_genome[get<0>(scored_genomes)[0]]]->write(record);
        }
        else if (get<0>(scored_genomes).size() == 0) {
            std::string record = ">";
            record.append(id);
            record.append("\n");
            record.append(seq);
            record.append("\n");
            fasta_writer["unmapped"]->write(record);
            _unmatched++;
        }
        else {
            // Ambiguous
            string header_tail = "|";
            for (auto const& _genome_id : get<0>(scored_genomes)) {
                header_tail.append(id_to_genome[_genome_id]);
                header_tail.append(";");
                stats[id_to_genome[_genome_id]]["ambig"]++;
            }
            header_tail.append(to_string(get<1>(scored_genomes)));
            std::string record = ">";
            record.append(id);
            record.append(header_tail);
            record.append("\n");
            record.append(seq);
            record.append("\n");
            fasta_writer["ambig"]->write(record);
        }

        total++;
        if (total == 5000) {
            cout << "processed: " << 5000 * ++chunks << "contigs" << endl;
            total = 0;
        }

    }

    for (auto const& stat : stats) {
        string genome_name = stat.first;
        auto _map = stat.second;
        int unique = _map["unique"];
        int ambig = _map["ambig"];
        cout << genome_name << ": unique(" << unique << ") " << "ambig(" << ambig << ")" << endl;
    }
    cout << "_unmatched: " << _unmatched << endl;

    for (auto f : fasta_writer)
        f.second->close();


    return 0;
}