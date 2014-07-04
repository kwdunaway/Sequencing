PATH=$PATH:/data/scratch/programs/tuxedo/bowtie-0.12.7/:/data/scratch/programs/tuxedo/samtools/

#SampleName=JLKD010 FASTQ_folder=Sample_JLKD010
gunzip -c Sample_JLKD010/*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > Sample_JLKD010/JLKD010.fq
split -l 80000000 Sample_JLKD010/JLKD010.fq Sample_JLKD010/JLKD010_
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD010/JLKD010_aa -o JLKD010_aa.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD010/JLKD010_ab -o JLKD010_ab.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD010/JLKD010_ac -o JLKD010_ac.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD010/JLKD010_ad -o JLKD010_ad.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD010/JLKD010_ae -o JLKD010_ae.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD010/JLKD010_af -o JLKD010_af.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD010/JLKD010_ag -o JLKD010_ag.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD010/JLKD010_ah -o JLKD010_ah.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD010/JLKD010_ai -o JLKD010_ai.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD010/JLKD010_aj -o JLKD010_aj.txt
cat JLKD010_a* > JLKD010_BSOUT_e80m2.txt
cat log_JLKD010_a* > log_JLKD010_BSOUT_e80m2.txt
perl /data/scratch/programs/perl_script/BS_Seeker_out2browserview_v4.pl JLKD010_BSOUT_e80m2.txt Sorted_JLKD010 JLKD010 hg18 n

#SampleName=JLKD011 FASTQ_folder=Sample_JLKD011
gunzip -c Sample_JLKD011/*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > Sample_JLKD011/JLKD011.fq
split -l 80000000 Sample_JLKD011/JLKD011.fq Sample_JLKD011/JLKD011_
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD011/JLKD011_aa -o JLKD011_aa.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD011/JLKD011_ab -o JLKD011_ab.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD011/JLKD011_ac -o JLKD011_ac.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD011/JLKD011_ad -o JLKD011_ad.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD011/JLKD011_ae -o JLKD011_ae.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD011/JLKD011_af -o JLKD011_af.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD011/JLKD011_ag -o JLKD011_ag.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD011/JLKD011_ah -o JLKD011_ah.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD011/JLKD011_ai -o JLKD011_ai.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD011/JLKD011_aj -o JLKD011_aj.txt
cat JLKD011_a* > JLKD011_BSOUT_e80m2.txt
cat log_JLKD011_a* > log_JLKD011_BSOUT_e80m2.txt
perl /data/scratch/programs/perl_script/BS_Seeker_out2browserview_v4.pl JLKD011_BSOUT_e80m2.txt Sorted_JLKD011 JLKD011 hg18 n

#SampleName=JLKD012 FASTQ_folder=Sample_JLKD012
gunzip -c Sample_JLKD012/*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > Sample_JLKD012/JLKD012.fq
split -l 80000000 Sample_JLKD012/JLKD012.fq Sample_JLKD012/JLKD012_
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD012/JLKD012_aa -o JLKD012_aa.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD012/JLKD012_ab -o JLKD012_ab.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD012/JLKD012_ac -o JLKD012_ac.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD012/JLKD012_ad -o JLKD012_ad.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD012/JLKD012_ae -o JLKD012_ae.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD012/JLKD012_af -o JLKD012_af.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD012/JLKD012_ag -o JLKD012_ag.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD012/JLKD012_ah -o JLKD012_ah.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD012/JLKD012_ai -o JLKD012_ai.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD012/JLKD012_aj -o JLKD012_aj.txt
cat JLKD012_a* > JLKD012_BSOUT_e80m2.txt
cat log_JLKD012_a* > log_JLKD012_BSOUT_e80m2.txt
perl /data/scratch/programs/perl_script/BS_Seeker_out2browserview_v4.pl JLKD012_BSOUT_e80m2.txt Sorted_JLKD012 JLKD012 hg18 n

#SampleName=JLKD013 FASTQ_folder=Sample_JLKD013
gunzip -c Sample_JLKD013/*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > Sample_JLKD013/JLKD013.fq
split -l 80000000 Sample_JLKD013/JLKD013.fq Sample_JLKD013/JLKD013_
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD013/JLKD013_aa -o JLKD013_aa.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD013/JLKD013_ab -o JLKD013_ab.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD013/JLKD013_ac -o JLKD013_ac.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD013/JLKD013_ad -o JLKD013_ad.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD013/JLKD013_ae -o JLKD013_ae.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD013/JLKD013_af -o JLKD013_af.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD013/JLKD013_ag -o JLKD013_ag.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD013/JLKD013_ah -o JLKD013_ah.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD013/JLKD013_ai -o JLKD013_ai.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD013/JLKD013_aj -o JLKD013_aj.txt
cat JLKD013_a* > JLKD013_BSOUT_e80m2.txt
cat log_JLKD013_a* > log_JLKD013_BSOUT_e80m2.txt
perl /data/scratch/programs/perl_script/BS_Seeker_out2browserview_v4.pl JLKD013_BSOUT_e80m2.txt Sorted_JLKD013 JLKD013 hg18 n

#SampleName=JLKD014 FASTQ_folder=Sample_JLKD014
gunzip -c Sample_JLKD014/*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > Sample_JLKD014/JLKD014.fq
split -l 80000000 Sample_JLKD014/JLKD014.fq Sample_JLKD014/JLKD014_
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD014/JLKD014_aa -o JLKD014_aa.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD014/JLKD014_ab -o JLKD014_ab.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD014/JLKD014_ac -o JLKD014_ac.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD014/JLKD014_ad -o JLKD014_ad.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD014/JLKD014_ae -o JLKD014_ae.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD014/JLKD014_af -o JLKD014_af.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD014/JLKD014_ag -o JLKD014_ag.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD014/JLKD014_ah -o JLKD014_ah.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD014/JLKD014_ai -o JLKD014_ai.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD014/JLKD014_aj -o JLKD014_aj.txt
cat JLKD014_a* > JLKD014_BSOUT_e80m2.txt
cat log_JLKD014_a* > log_JLKD014_BSOUT_e80m2.txt
perl /data/scratch/programs/perl_script/BS_Seeker_out2browserview_v4.pl JLKD014_BSOUT_e80m2.txt Sorted_JLKD014 JLKD014 hg18 n

#SampleName=JLKD015 FASTQ_folder=Sample_JLKD015
gunzip -c Sample_JLKD015/*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > Sample_JLKD015/JLKD015.fq
split -l 80000000 Sample_JLKD015/JLKD015.fq Sample_JLKD015/JLKD015_
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD015/JLKD015_aa -o JLKD015_aa.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD015/JLKD015_ab -o JLKD015_ab.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD015/JLKD015_ac -o JLKD015_ac.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD015/JLKD015_ad -o JLKD015_ad.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD015/JLKD015_ae -o JLKD015_ae.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD015/JLKD015_af -o JLKD015_af.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD015/JLKD015_ag -o JLKD015_ag.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD015/JLKD015_ah -o JLKD015_ah.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD015/JLKD015_ai -o JLKD015_ai.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD015/JLKD015_aj -o JLKD015_aj.txt
cat JLKD015_a* > JLKD015_BSOUT_e80m2.txt
cat log_JLKD015_a* > log_JLKD015_BSOUT_e80m2.txt
perl /data/scratch/programs/perl_script/BS_Seeker_out2browserview_v4.pl JLKD015_BSOUT_e80m2.txt Sorted_JLKD015 JLKD015 hg18 n

#SampleName=JLKD016 FASTQ_folder=Sample_JLKD016
gunzip -c Sample_JLKD016/*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > Sample_JLKD016/JLKD016.fq
split -l 80000000 Sample_JLKD016/JLKD016.fq Sample_JLKD016/JLKD016_
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD016/JLKD016_aa -o JLKD016_aa.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD016/JLKD016_ab -o JLKD016_ab.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD016/JLKD016_ac -o JLKD016_ac.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD016/JLKD016_ad -o JLKD016_ad.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD016/JLKD016_ae -o JLKD016_ae.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD016/JLKD016_af -o JLKD016_af.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD016/JLKD016_ag -o JLKD016_ag.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD016/JLKD016_ah -o JLKD016_ah.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD016/JLKD016_ai -o JLKD016_ai.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD016/JLKD016_aj -o JLKD016_aj.txt
cat JLKD016_a* > JLKD016_BSOUT_e80m2.txt
cat log_JLKD016_a* > log_JLKD016_BSOUT_e80m2.txt
perl /data/scratch/programs/perl_script/BS_Seeker_out2browserview_v4.pl JLKD016_BSOUT_e80m2.txt Sorted_JLKD016 JLKD016 hg18 n

#SampleName=JLKD017 FASTQ_folder=Sample_JLKD017
gunzip -c Sample_JLKD017/*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > Sample_JLKD017/JLKD017.fq
split -l 80000000 Sample_JLKD017/JLKD017.fq Sample_JLKD017/JLKD017_
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD017/JLKD017_aa -o JLKD017_aa.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD017/JLKD017_ab -o JLKD017_ab.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD017/JLKD017_ac -o JLKD017_ac.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD017/JLKD017_ad -o JLKD017_ad.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD017/JLKD017_ae -o JLKD017_ae.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD017/JLKD017_af -o JLKD017_af.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD017/JLKD017_ag -o JLKD017_ag.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD017/JLKD017_ah -o JLKD017_ah.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD017/JLKD017_ai -o JLKD017_ai.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD017/JLKD017_aj -o JLKD017_aj.txt
cat JLKD017_a* > JLKD017_BSOUT_e80m2.txt
cat log_JLKD017_a* > log_JLKD017_BSOUT_e80m2.txt
perl /data/scratch/programs/perl_script/BS_Seeker_out2browserview_v4.pl JLKD017_BSOUT_e80m2.txt Sorted_JLKD017 JLKD017 hg18 n

#SampleName=JLKD018 FASTQ_folder=Sample_JLKD018
gunzip -c Sample_JLKD018/*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > Sample_JLKD018/JLKD018.fq
split -l 80000000 Sample_JLKD018/JLKD018.fq Sample_JLKD018/JLKD018_
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD018/JLKD018_aa -o JLKD018_aa.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD018/JLKD018_ab -o JLKD018_ab.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD018/JLKD018_ac -o JLKD018_ac.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD018/JLKD018_ad -o JLKD018_ad.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD018/JLKD018_ae -o JLKD018_ae.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD018/JLKD018_af -o JLKD018_af.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD018/JLKD018_ag -o JLKD018_ag.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD018/JLKD018_ah -o JLKD018_ah.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD018/JLKD018_ai -o JLKD018_ai.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD018/JLKD018_aj -o JLKD018_aj.txt
cat JLKD018_a* > JLKD018_BSOUT_e80m2.txt
cat log_JLKD018_a* > log_JLKD018_BSOUT_e80m2.txt
perl /data/scratch/programs/perl_script/BS_Seeker_out2browserview_v4.pl JLKD018_BSOUT_e80m2.txt Sorted_JLKD018 JLKD018 hg18 n

#SampleName=JLKD019 FASTQ_folder=Sample_JLKD019
gunzip -c Sample_JLKD019/*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > Sample_JLKD019/JLKD019.fq
split -l 80000000 Sample_JLKD019/JLKD019.fq Sample_JLKD019/JLKD019_
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD019/JLKD019_aa -o JLKD019_aa.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD019/JLKD019_ab -o JLKD019_ab.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD019/JLKD019_ac -o JLKD019_ac.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD019/JLKD019_ad -o JLKD019_ad.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD019/JLKD019_ae -o JLKD019_ae.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD019/JLKD019_af -o JLKD019_af.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD019/JLKD019_ag -o JLKD019_ag.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD019/JLKD019_ah -o JLKD019_ah.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD019/JLKD019_ai -o JLKD019_ai.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD019/JLKD019_aj -o JLKD019_aj.txt
cat JLKD019_a* > JLKD019_BSOUT_e80m2.txt
cat log_JLKD019_a* > log_JLKD019_BSOUT_e80m2.txt
perl /data/scratch/programs/perl_script/BS_Seeker_out2browserview_v4.pl JLKD019_BSOUT_e80m2.txt Sorted_JLKD019 JLKD019 hg18 n

#SampleName=JLKD020 FASTQ_folder=Sample_JLKD020
gunzip -c Sample_JLKD020/*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > Sample_JLKD020/JLKD020.fq
split -l 80000000 Sample_JLKD020/JLKD020.fq Sample_JLKD020/JLKD020_
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD020/JLKD020_aa -o JLKD020_aa.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD020/JLKD020_ab -o JLKD020_ab.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD020/JLKD020_ac -o JLKD020_ac.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD020/JLKD020_ad -o JLKD020_ad.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD020/JLKD020_ae -o JLKD020_ae.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD020/JLKD020_af -o JLKD020_af.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD020/JLKD020_ag -o JLKD020_ag.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD020/JLKD020_ah -o JLKD020_ah.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD020/JLKD020_ai -o JLKD020_ai.txt
python /data/scratch/programs/tuxedo/BS_Seeker/BS_Seeker.py -t N -e 80 -m 2 -a /data/scratch/programs/tuxedo/BS_Seeker/adapter2.txt -d /data/scratch/genomes/hg18/reference_genome/ -p /data/scratch/programs/tuxedo/bowtie-0.12.7/ -i Sample_JLKD020/JLKD020_aj -o JLKD020_aj.txt
cat JLKD020_a* > JLKD020_BSOUT_e80m2.txt
cat log_JLKD020_a* > log_JLKD020_BSOUT_e80m2.txt
perl /data/scratch/programs/perl_script/BS_Seeker_out2browserview_v4.pl JLKD020_BSOUT_e80m2.txt Sorted_JLKD020 JLKD020 hg18 n

