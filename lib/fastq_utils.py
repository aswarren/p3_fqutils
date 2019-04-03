#!/usr/bin/env python

import os, sys
import subprocess
import multiprocessing
import tarfile, json
import requests

#take genome data structure and make directory names.
def make_directory_names(genome):
    if genome["dir"].endswith('/'):
        genome["dir"]=genome["dir"][:-1]
    genome["dir"]=os.path.abspath(genome["dir"])
    genome["output"]=os.path.join(output_dir,os.path.basename(genome["dir"]))

#hisat2 has problems with spaces in filenames
#prevent spaces in filenames. if one exists link the file to a no-space version.
def link_space(file_path):
    result=file_path
    name=os.path.splitext(os.path.basename(file_path))[0]
    if " " in name:
        clean_name=name.replace(" ","")
        result= file_path.replace(name,clean_name)
        if not os.path.exists(result):
            subprocess.check_call(["ln","-s",file_path,result])
    return result


def run_fastqc(read_list, output_dir, job_data):
    rcount=0
    for r in read_list:
        cur_cleanup=[]
        rcount+=1
        fastqc_cmd=["fastqc","--outdir",output_dir]
        cur_cmd=list(cmd)
        if "read2" in r:
            fastqc_cmd+=[r["read1"],r["read2"]]
        else:
            fastqc_cmd+=[r["read1"]]
        print " ".join(fastqc_cmd)
        subprocess.check_call(fastqc_cmd)

def run_alignment(genome_list, read_list, parameters, output_dir, job_data): 
    #modifies condition_dict sub replicates to include 'bowtie' dict recording output files
    for genome in genome_list:
        genome_link = genome["genome_link"]
        final_cleanup=[]
        if "hisat_index" in genome and genome["hisat_index"]:
            archive = tarfile.open(genome["hisat_index"])
            indices= [os.path.join(output_dir,os.path.basename(x)) for x in archive.getnames()]
            final_cleanup+=indices
            #archive.extractall(path=output_dir)
            archive.close()
            subprocess.check_call(["tar","-xvf", genome["hisat_index"], "-C", output_dir])
            index_prefix = os.path.join(output_dir, os.path.basename(genome["hisat_index"]).replace(".ht2.tar","")) #somewhat fragile convention. tar prefix is underlying index prefix
            cmd=["hisat2","--dta-cufflinks", "-x", index_prefix] 
            thread_count= parameters.get("hisat2",{}).get("-p",0)
        else:
            subprocess.check_call(["bowtie2-build", genome_link, genome_link])
            #cmd=["hisat2","--dta-cufflinks", "-x", genome_link, "--no-spliced-alignment"] 
            cmd=["bowtie2", "-x", genome_link]
            thread_count= parameters.get("bowtie2",{}).get("-p",0)
        if thread_count == 0:
            thread_count=2 #multiprocessing.cpu_count()
        cmd+=["-p",str(thread_count)]
        target_dir=genome['output']
        for r in read_list:
            rcount=0
            cur_cleanup=[]
            rcount+=1
            samstat_cmd=["samstat"]
            cur_cmd=list(cmd)
            if "read2" in r:
                cur_cmd+=["-1",link_space(r["read1"])," -2",link_space(r["read2"])]
                name1=os.path.splitext(os.path.basename(r["read1"]))[0].replace(" ","")
                name2=os.path.splitext(os.path.basename(r["read2"]))[0].replace(" ","")
                sam_file=os.path.join(target_dir,name1+"_"+name2+".sam")
            else:
                cur_cmd+=[" -U",link_space(r["read1"])]
                name1=os.path.splitext(os.path.basename(r["read1"]))[0].replace(" ","")
                sam_file=os.path.join(target_dir,name1+".sam")
            cur_cleanup.append(sam_file)
            bam_file=sam_file[:-4]+".bam"
            samstat_cmd.append(bam_file)
            r[genome["genome"]]={}
            r[genome["genome"]]["bam"]=bam_file
            cur_cmd+=["-S",sam_file]
            if os.path.exists(bam_file):
                sys.stderr.write(bam_file+" alignments file already exists. skipping\n")
            else:
                print cur_cmd
                subprocess.check_call(cur_cmd) #call bowtie2
            if not os.path.exists(bam_file):
                subprocess.check_call("samtools view -Su "+sam_file+" | samtools sort -o - - > "+bam_file, shell=True)#convert to bam
                subprocess.check_call("samtools index "+bam_file, shell=True)
                #subprocess.check_call('samtools view -S -b %s > %s' % (sam_file, bam_file+".tmp"), shell=True)
                #subprocess.check_call('samtools sort %s %s' % (bam_file+".tmp", bam_file), shell=True)
            print " ".join(samstat_cmd)
            subprocess.check_call(samstat_cmd)

            for garbage in cur_cleanup:
                subprocess.call(["rm", garbage])
        for garbage in final_cleanup:
            subprocess.call(["rm", garbage])

def get_genome(parameters):
    target_file = os.path.join(parameters["output_path"],parameters["gid"]+".fna")
    if not os.path.exists(target_file):
        genome_url= "data_url/genome_sequence/?eq(genome_id,gid)&limit(25000)".replace("data_url",parameters["data_url"]).replace("gid",parameters["gid"])
        print genome_url
        headers = {"accept":"application/sralign+dna+fasta"}
        #print "switch THE HEADER BACK!"
        #headers = {'Content-Type': 'application/x-www-form-urlencoded; charset=utf-8'}
        req = requests.Request('GET', genome_url, headers=headers)
        prepared = req.prepare()
        #pretty_print_POST(prepared)
        s = requests.Session()
        response=s.send(prepared)
        handle = open(target_file, 'w')
        if not response.ok:
            sys.stderr.write("API not responding. Please try again later.\n")
            sys.exit(2)
        else:
            for block in response.iter_content(1024):
                handle.write(block)
    return target_file


def setup(job_data, output_dir, tool_params):
    genome_ids=[job_data.get("reference_genome_id",None)]
    genome_list=[]
    for gid in genome_ids:
        genome={}
        job_data["gid"]=gid #cheat
        genome["genome_link"]=get_genome(job_data)
        genome={"gid"}=gid
        genome["output"=output_dir
        genome_list.append(genome)

    read_list = []
    for read in job_data.get("paired_end_libs",[])+job_data.get("single_end_libs",[])+job_data.get("srr_libs",[]):
        if "read" in read:
            read["read1"] = read.pop("read")
        read_list.append(read)

    rcount=0
    for r in read_list:
        target_dir=output_dir
        #subprocess.call(["mkdir","-p",target_dir])
        rcount+=1
        if "srr_accession" in r:
            srr_id = r["srr_accession"] 
            meta_file = os.path.join(target_dir,srr_id+"_meta.txt")
            subprocess.check_call(["p3-sra","--gzip","--out",target_dir,"--metadata-file", meta_file, "--id",srr_id])
            with open(meta_file) as f:
                job_meta = json.load(f)
                files = job_meta[0].get("files",[])
                if len(files) > 0:
                    for i,f in enumerate(files):
                        r["read"+str(i+1)]=os.path.join(target_dir, f)
    recipe = job_data.get("recipe",[])
    return genome_list, read_list, recipe






def run_fq_util(job_data, output_dir, tool_params):
    #arguments:
    #list of genomes [{"genome":somefile,"annotation":somefile}]
    #dictionary of library dictionaries structured as {libraryname:{library:libraryname, replicates:[{read1:read1file, read2:read2file}]}}
    #parametrs_file is json parameters list keyed as bowtie, cufflinks, cuffdiff.
    output_dir=os.path.abspath(output_dir)
    subprocess.call(["mkdir","-p",output_dir])

    genome_list, read_list, recipe=setup(job_data, output_dir, tool_params)
    for step in recipe:
        if step == "TRIM":
            trimmed_reads = run_trim(read_list, output_dir, job_data)
            read_list = trimmed_reads
        if step == "FASTQC":
            run_fastqc(read_list, output_dir, job_data)
        if step == "ALIGN":
            run_alignment(genome_list, read_list, parameters, output_dir, job_data)


