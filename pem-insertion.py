import logging
import os
import argparse
import sys
import pysam
import pandas as pd

def check_dir(dir):
    if not os.path.exists(dir):
        os.system('mkdir %s' % dir)

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq_rc = "".join(complement.get(base, base) for base in reversed(seq))
    return seq_rc

# def exec plus log function
def exec_log(cmdline):
    logging.info(cmdline)
    os.system(cmdline)

# prepare args and make output file
def parse_args():
    parser = argparse.ArgumentParser(description='This pipeline performs n bp insertion analysis.')
    parser.add_argument("-r1", type = str, required = True, help = "read 1 file" )
    parser.add_argument("-r2", type = str, required = True, help = "read 2 file" )
    parser.add_argument("-b", dest = "barcode", type = str, required = True, help = "library barcode" )
    parser.add_argument("-p", dest = "primer", type = str, required = True, help = "library primer" )
    parser.add_argument("-a", dest = "adapter", type = str, required = True, help = "library adapter" )
    parser.add_argument("-n", dest = "name", type = str, required = True, help = "library name" )
    parser.add_argument("-m", dest = "metafile", type = str, required = True, help = "metadata file documentation describe the translocation type")
    parser.add_argument("-o", dest = "outdir", type = str, required = True, help = "output directory" )
    parser.add_argument("-num", dest = "number", type =int, required = True, help = "detect n bp insertion")
    args = parser.parse_args()
    check_dir(args.outdir)

    return args

# create logfile
def create_logfile(args):
    logdir = args.outdir + '/pipeline.log'
    # init logging
    logging.basicConfig( level = 10,
    format = '%(levelname)-5s @ %(asctime)s: %(message)s ',
    datefmt = '%a, %d %b %Y %H:%M:%S',
    filename = logdir,
    filemode = 'a' )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(levelname)-5s @ %(asctime)s: %(message)s ','%a, %d %b %Y %H:%M:%S')
    #formatter.formatTime('%a, %d %b %Y %H:%M:%S')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    if not logging.getLogger('').handlers:
        logging.getLogger('').addHandler(console)
    else:
        logging.getLogger('')

# dedup the r1 and r2 reads based on the rmb sequence in the r2 reads
def rmb_dedup(fdedup,args):
    logging.info('Deduping fastq file by RMB ......')

    # make adapter.fa file
    adapter = args.adapter.strip()
    fadapter = open('{}/adapter.fa'.format(fdedup),"w")
    fadapter.write(">adapter"+"\n"+adapter+"\n")
    fadapter.close()

    # load fasta file

    os.system('gzip -c -d %s > %s/raw.r1.fq' % (args.r1, fdedup))
    os.system('gzip -c -d %s > %s/raw.r2.fq' % (args.r2, fdedup))

    # use bwa mem to map reads against adapter.fa file

    os.system("samtools faidx {0}/adapter.fa".format(fdedup))
    os.system("bwa index -p {0}/adapter {0}/adapter.fa 1>{0}/build_index.o 2>{0}/build_index.e".format(fdedup))
    os.system("bwa mem -t 4 {0}/adapter -k 10 -L 0 -T 10 {0}/raw.r2.fq > {0}/{1}.sam 2>{0}/bwa_align_adapter.log".format(fdedup,args.name))
    os.system("samtools view -S -b -h {0}/{1}.sam > {0}/{1}.bam".format(fdedup,args.name))
    os.system("samtools sort -@ 4 {0}/{1}.bam > {0}/{1}.sort.bam".format(fdedup,args.name))
    os.system("samtools index {0}/{1}.sort.bam".format(fdedup,args.name))

    # dedup based on the bam file, get the uniq reads ID
    bam_file = pysam.AlignmentFile("{0}/{1}.sort.bam".format(fdedup,args.name), "rb")

    barcode_list = open("{0}/{1}_barcode_list.txt".format(fdedup,args.name), "w")

    for read in bam_file.fetch():
        if read.query_alignment_length >= (16 - 2):
            bar_start = read.query_alignment_end + 1
            bar_end = read.query_length
            bar_length = bar_end - bar_start + 1

            if bar_length <= 17 and bar_length >= (17 - 2):
                barcode_seq = read.query_sequence[(bar_start - 1): bar_end]
                barcode_list.write(str(read.query_name) + "\t" + str(barcode_seq) + "\n")

    data = pd.read_csv("{0}/{1}_barcode_list.txt".format(fdedup,args.name), sep='\t', names=["Qname", "Barcode"])

    # get barcode length
    length = data['Barcode'].astype('str')
    length = data['Barcode'].str.len()

    data['Length'] = length

    # get barcode frequency
    data['Freq'] = data.groupby('Barcode')['Barcode'].transform('count')
    data = data.sort_values(by=['Freq', 'Barcode', 'Length'], ascending=False)
    data.to_csv("{0}/{1}_barcode_sort.txt".format(fdedup,args.name), header=True, sep='\t', index=False,
                columns=[u'Qname', u'Barcode', u'Freq', u'Length'])

    # also generate duplicated barcode
    datb = data[data.duplicated(['Barcode'])]
    datb = datb.reset_index(drop=True)

    # unique barcode
    data = data.drop_duplicates(['Barcode'])  # rough dedup
    data = data.reset_index(drop=True)  # reset index

    data.to_csv("{0}/{1}_barcode_uniq.txt".format(fdedup,args.name), header=True, sep='\t', index=False,
                columns=[u'Qname', u'Barcode', u'Freq', u'Length'])
    datb.to_csv("{0}/{1}_barcode_dup.txt".format(fdedup,args.name), header=True, sep='\t', index=False,
                columns=[u'Qname', u'Barcode', u'Freq', u'Length'])

    os.system("less {0}/{1}_barcode_uniq.txt|cut -f1|sort > {0}/{1}_uniq_ID.list ".format(fdedup,args.name))

    logging.info('Get the unique reads ID.')
    # uniq raw r1 r2 fastq file
    exec_log('seqtk subseq {0}/raw.r1.fq {0}/{1}_uniq_ID.list > {0}/{1}_R1.uniq.fq'.format(fdedup,args.name))
    exec_log('seqtk subseq {0}/raw.r2.fq {0}/{1}_uniq_ID.list > {0}/{1}_R2.uniq.fq'.format(fdedup,args.name))
    logging.info('Get the deduped r1 r2 fastq file.')
    # sort the directory
    os.system('mkdir {0}/index'.format(fdedup))
    os.system('mv {0}/adapter* {0}/index'.format(fdedup))
    os.system('rm -rf {0}/raw.r1.fq {0}/raw.r2.fq {0}/{1}.sam {0}/{1}.bam'.format(fdedup,args.name))

def stitch(fdedup,ffastq,args):
    logging.info("Joining the uniq r1 and r2 reads ......")
    # -o prefix -d directory
    exec_log('flash {0}/{1}_R1.uniq.fq {0}/{1}_R2.uniq.fq -o {1} -d {2}/flash_out > {2}/flash.log'.format(fdedup,args.name,ffastq))
    fstitch=ffastq+"/flash_out/"+args.name
    f1 = fstitch + '.extendedFrags.fastq'
    f2 = fstitch + '.notCombined_1.fastq'
    merge_f = fstitch + '.merge.fastq'
    cmd = "cat {} {} > {}".format(f1, f2, merge_f)
    os.system(cmd)
    cmd = "gzip {}".format(fstitch + '.merge.fastq')
    os.system(cmd)
    logging.info('Get merged fastq file! ')

def demutiplex_trim_barcode(ffastq,args):
    barcode = args.barcode.strip()
    primer=args.primer.strip()
    de_barcode=(barcode+primer)[:12]
    adapter = args.adapter.strip()
    adapter_rc = reverse_complement(adapter)

    # prepare de multiplex barcode.txt
    fbar=open('{}/barcode.txt'.format(ffastq),"w")
    fbar.write('{}\t{}\n'.format(args.name,de_barcode))
    fbar.close()

    # de multiplex
    exec_log('fastq-multx -m 0 -x -b -d 0 -B {0}/barcode.txt {0}/flash_out/{1}.merge.fastq.gz -o {0}/%_demultiplex_merge.fq >{0}/fastq-multx.log 2>&1'.format(ffastq, args.name))
    os.system('mkdir {0}/unmatched'.format(ffastq))
    os.system('mv {0}/unmatched_* {0}/unmatched'.format(ffastq))

    # trim merge
    exec_log('cutadapt -g %s -n 1 -m 50 -o %s/%s_demultiplex_trim_tmp_merge.fq %s/%s_demultiplex_merge.fq >%s/cutadaper1.log 2>&1' % (
    barcode, ffastq, args.name, ffastq, args.name,ffastq))
    exec_log('cutadapt -a %s -n 1 -m 50 -o %s/%s_demultiplex_trim_merge.fq %s/%s_demultiplex_trim_tmp_merge.fq >%s/cutadaper2.log 2>&1' % (
    adapter, ffastq, args.name, ffastq, args.name,ffastq))
# cp meta_file to outdir
def check_metafile(args):
    try:
        os.system("perl -p -e 's/\\r/\\n/g' %s > %s/metadata.txt" % (args.metafile, args.outdir))
    except:
        logging.error('Problem converting line feeds to unix format')
        sys.exit(-1)

# load information in meta_file into a directory
def load_metafile(args):
    # Parse metafile
    metapath=args.outdir+'/metadata.txt'
    si={}
    for line in open(metapath):
        if line.startswith('translocation'): continue
        l = line.split(',')
        try:
            samplename = '%s' % (l[0])
            si[samplename] = [l[1].upper().strip(), l[2].strip(), l[3].strip()]
        except:
            logging.error('Cannot find primer and adapter sequences in meta file')
            sys.exit(-1)
    return si

def make_result_dic(n):
    import itertools
    all_type=[]
    for i in itertools.product("ATGC", repeat=n):
        all_type.append("".join(i))
    result_dic={i:0 for i in all_type}
    return (all_type, result_dic)

def map_against_translocation_reference(si,sample,fsample,ffastq,args):

    # make translocation.fa file
    translocation_name=sample
    translocation_fasta=si[sample][0].strip()
    if not os.path.exists('{0}/{1}.fa'.format(fsample,translocation_name)):
        logging.info("making own reference sequence")
        ftfasta = open('{0}/{1}.fa'.format(fsample,translocation_name),"w")
        ftfasta.write(">"+translocation_name+"\n"+translocation_fasta+"\n")
        ftfasta.close()

    # bwa-mem map merged reads against translocation.fa file
    if not os.path.exists('{0}/{1}.sort.bam'.format(fsample,translocation_name)):
        logging.info('map reads against sequence of own reference sequence')
        os.system("samtools faidx {0}/{1}.fa > {0}/samtools_index.log 2>&1 ".format(fsample,translocation_name))
        os.system("bwa index -p {0}/{1} {0}/{1}.fa 1>{0}/build_index.o 2>{0}/build_index.e".format(fsample,translocation_name))
        exec_log("bwa mem -t 4 {0}/{1} -k 10 -L 0 -T 10 {2}/{3}_demultiplex_trim_merge.fq > {0}/{1}.sam 2>{0}/bwa_align_adapter.log".format(fsample,translocation_name,ffastq,args.name))
        os.system("samtools view -S -b -h {0}/{1}.sam > {0}/{1}.bam".format(fsample,translocation_name))
        os.system("samtools sort -@ 4 {0}/{1}.bam > {0}/{1}.sort.bam".format(fsample,translocation_name))
        os.system("samtools index {0}/{1}.sort.bam".format(fsample,translocation_name))
        os.system("rm -rf {0}/{1}.sam {0}/{1}.bam".format(fsample, translocation_name))

def define_insersion(fsample,sample,si,n,all_type,result_dic,args):
    # identify the sequence fit
    translocation_name = sample
    bait_length = int(si[sample][1].strip())
    unique_length = int(si[sample][2].strip())

    if not os.path.exists("{0}/{1}_{2}bp_insertion_pic.txt".format(fsample,translocation_name,n)):
        # get bait_end and prey start sequence
        for line in open("{0}/{1}.fa".format(fsample, translocation_name), "r"):
            if not line.startswith(">"):
                # eg.bait length=72 ,unique_length =20, bait end[52:72], prey start[72:92]
                bait_end = line[bait_length - unique_length: bait_length].strip()
                prey_start = line[bait_length:bait_length + unique_length].strip()
                logging.info('bait_end: {}; prey_start:{}'.format(bait_end, prey_start))

        bam_file = pysam.AlignmentFile("{0}/{1}.sort.bam".format(fsample, translocation_name), "rb")
        output_file = open("{0}/{1}_{2}bp_insertion.txt".format(fsample, translocation_name, n), "w")

        # based on bait length and unique length information get the bait_end and prey_start sequence
        for read in bam_file.fetch():
            for insert in all_type:

                unique_query = bait_end.strip() + insert + prey_start.strip()
                test = unique_query in read.query_sequence
                if test == True:

                    list = read.cigartuples
                    list_dic = {i: 0 for i in range(0, 10)}

                    for i in list:
                        list_dic[i[0]] += i[1]

                    if list_dic[4] < 10 and list_dic[2] < 10:
                        result_dic[insert] += 1
                        output_file.write(str(read.query_name) + "\t" + str(read.cigarstring) + "\t"
                                          + str(read.query_sequence)[bait_length:bait_length + n] + "\n")

        output_file.close()

        of = open("{0}/{1}_{2}bp_insertion_count.txt".format(fsample, translocation_name, n), "w")

        for i in result_dic:
            of.write(i + "\t" + str(result_dic[i]) + "\n")
        of.close()

        of1 = open("{0}/{1}_{2}bp_insertion_pic.txt".format(fsample, translocation_name, n), "w")
        for line in open("{0}/{1}_{2}bp_insertion_count.txt".format(fsample, translocation_name, n), "r"):
            l = line.strip().split("\t")
            if int(l[1]) > 0:
                of1.write(line.strip() + "\n")
        of1.close()

        check_dir('{0}/result'.format(args.outdir))
        os.system('cp {0}/{1}_{2}bp_insertion_count.txt {3}/result'.format(fsample, translocation_name, n, args.outdir))
        logging.info("{0}/{1}_{2}bp_insertion_count.txt get, nice job, see you~".format(fsample, sample, n))

def main():
    args = parse_args()
    create_logfile(args)

    logging.info('Welcome to inbp pipeline!')
    logging.info('Paramaters: ' + ' '.join(sys.argv))
# dedup raw fastq
    fdedup = '{}/dedup'.format(args.outdir)
    # check if already done
    if os.path.exists('{}'.format(fdedup)):
        logging.info('Already deduped the raw fastq file!')
    else:
        os.system('mkdir {}'.format(fdedup))
        rmb_dedup(fdedup, args)

# prepare fastq for map
    ffastq = '{}/fastq'.format(args.outdir)
    # check if already done
    if os.path.exists('{}'.format(ffastq)):
        logging.info('Already prepared the unique fastq file!')
    else:
        os.system('mkdir {}'.format(ffastq))
        stitch(fdedup, ffastq, args)
        demutiplex_trim_barcode(ffastq, args)

    check_metafile(args)
    si=load_metafile(args)

    # make identify dir to put all the translocation result file
    fidentify= '{}/identify'.format(args.outdir)
    check_dir(fidentify)

    # make every different translocation bam file
    for sample in si:
        fsample = '{}/{}'.format(fidentify, sample)
        check_dir(fsample)
        map_against_translocation_reference(si, sample, fsample, ffastq, args)

        for n in range(1,args.number+1):
            all_type = make_result_dic(n)[0]
            result_dic = make_result_dic(n)[1]
            define_insersion(fsample, sample, si, n, all_type, result_dic, args)

        os.system("cat {3}/result/{1}_*bp_insertion_count.txt > {3}/{1}_insertion_count.txt".format(fsample, sample, args.number, args.outdir))

if __name__ == '__main__':
    main()
main()
