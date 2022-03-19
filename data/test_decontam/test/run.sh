mkdir sketches
cd sketches
kSpider sketch --fastx ../cont1.fasta -k 21
kSpider sketch --fastx ../cont2.fasta -k 21
cd ..
reads=10_transcript.fasta

/home/mabuelanin/dib-dev/refDecontam/build/decontam_by_sketches $reads sketches

