# import os
# import psutil
# import shlex
# import time
#
# script_time = time.strftime( "%d/%m/%Y %H:%M:%S", time.localtime( time.time( ) ) )
# print("START %s" % script_time)
#
#
# def pid_status(process_pid):
#     status = ""
#     if psutil.pid_exists( process_pid ):
#         status = psutil.Process( process_pid ).status( )
#     else:
#         status = "finished"
#     return status
#
# def monitoring(command: str, info_store: dict):
#     cmd = shlex.split(command)
#     p = psutil.Popen(cmd)
#     process_id = p.pid
#     while pid_status(p.pid) not in ["finished", "zombie"]:
#         name = "%s####%s###%s" % (str(process_id), " ".join(p.cmdline()), str(p.create_time()))
#         info_store.setdefault(name, {})
#         info_store[name].setdefault("memory", set())
#         info_store[name]["memory"].add(str(p.memory_info()[0]))
#         info_store[name].setdefault("threads", set())
#         info_store[name]["threads"].add(str(p.num_threads()))
#         for children in p.children(recursive=True):
#             children_pid = children.pid
#             if psutil.pid_exists(children_pid):
#                 try:
#                     name = "%s####%s###%s" % (
#                     str(children_pid), " ".join(children.cmdline()), str(children.create_time()))
#                     info_store.setdefault(name, {})
#                     info_store[name].setdefault("memory", set())
#                     info_store[name]["memory"].add(str(children.memory_info()[0]))
#                     info_store[name].setdefault("threads", set())
#                     info_store[name]["threads"].add(str(children.num_threads()))
#                 except:
#                     pass
#     return info_store
import os.path
import cmdbench
import pandas as pd


if __name__ == '__main__':
    # outPath = "/home3/gdefazio/kMetaShot_benchmark/cami2/time_memo_assess/air_sr_s12"
    outPath = "/home3/gdefazio/kMetaShot_benchmark/cami2/time_memo_assess/air_pb_s4_ts"
    # samplePath="/home3/gdefazio/kMetaShot_benchmark/cami2/Airways/short_read/2017.12.04_18.56.22_sample_12/contigs"
    samplePath = "/home3/gdefazio/kMetaShot_benchmark/cami2/Airways/pacbio/2018.01.23_14.08.31_sample_4/contigs"
    kmetashotREF="/home3/gdefazio/kMetaShot_Bacteria_Archaea_ref/clndMashOne4strain/compressed_clndMashOne4strain_bacteria_archaea.h5"
    kmetashot_command= """
    conda run -n gius \
    /home/gdefazio/miniconda2/envs/gius/bin/kMetaShot_classifier_NV.py \
    -b %s/bins/ \
    -r %s \
    -p 10 \
    -a 0.1 \
    -o %s/kMetaShot_classif""" % (samplePath, kmetashotREF, outPath)

    camitaxREF = "/home3/gdefazio/kMetaShot_benchmark/camitax_reference/"

    camitax_command = """
    conda run -n camitax /home/gdefazio/miniconda2/envs/camitax/bin/nextflow \
    -log ./camitax.log \
    -bg run CAMI-challenge/CAMITAX \
    -with-trace \
    -with-timeline timeline.txt \
    --db %s \
    --i %s/bins  \
    --x fa""" % (camitaxREF, samplePath)

    gtdbtk_command = """
    conda run -n gtdbtk 
    /home/gdefazio/miniconda2/envs/gtdbtk/bin/gtdbtk \
    classify_wf  --genome_dir %s/bins \
    --out_dir %s/gtdbtk_classif --cpus 10 -x fa""" % (samplePath, outPath)

    # print("\n\nKMETASHOT\n")
    # kmetashot_bench = cmdbench.benchmark_command(kmetashot_command, iterations_num = 1)
    # kmetashot_bench_ref = kmetashot_bench.get_first_iteration()
    # print(kmetashot_bench_ref)
    # # pd.DataFrame(kmetashot_bench_ref['time_series']).to_csv(os.path.join(outPath,'kms_air_sr_s12_ts.csv'))
    # pd.DataFrame(kmetashot_bench_ref['time_series']).to_csv(os.path.join(outPath,'kms_air_pb_s4_ts.csv'))

    print("\n\nCAMITAX\n")
    camitax_bench = cmdbench.benchmark_command(camitax_command, iterations_num=1)
    camitax_bench_ref = camitax_bench.get_first_iteration()
    print(camitax_bench_ref)
    pd.DataFrame(camitax_bench_ref['time_series']).to_csv(os.path.join(outPath, 'cmx_air_sr_s12_ts.csv'))
    # pd.DataFrame(camitax_bench_ref['time_series']).to_csv(os.path.join(outPath, 'cmx_air_pb_s4_ts.csv'))

    # print("\n\nGTDBTK\n")
    # gtdbtk_bench = cmdbench.benchmark_command(gtdbtk_command, iterations_num=1)
    # gtdbtk_bench_ref = gtdbtk_bench.get_first_iteration()
    # print(gtdbtk_bench_ref)
    # pd.DataFrame(gtdbtk_bench_ref['time_series']).to_csv(os.path.join(outPath, 'gtk_air_sr_s12_ts.csv'))
    # pd.DataFrame(gtdbtk_bench_ref['time_series']).to_csv(os.path.join(outPath, 'gtk_air_pb_s4_ts.csv'))



    # with open("technical_measurement.tsv", "w") as tmp:
    #     tmp.write("process\tmemory\tthreads\n")
    #     for process in process_info_consumtion.keys():
    #         tmp.write("%s\t%s\t%s\n" % (
    #             process, ";".join(list(process_info_consumtion[process]["memory"])),
    #             ";".join(list(process_info_consumtion[process]["threads"]))))


# # kMETASHOT
# cmd = shlex.split( "python testing_script_folder/MetaShot_Master_script.py -p parameters_file.txt -m read_list" )
# p = psutil.Popen( cmd )
# process_id = p.pid
# while pid_status( p.pid ) not in ["finished", "zombie"]:
#     name = "%s####%s###%s" % (str( process_id ), " ".join( p.cmdline( ) ), str( p.create_time( ) ))
#     process_info_consumtion.setdefault( name, {} )
#     process_info_consumtion[name].setdefault( "memory", set( ) )
#     process_info_consumtion[name]["memory"].add( str( p.memory_info( )[0] ) )
#     process_info_consumtion[name].setdefault( "threads", set( ) )
#     process_info_consumtion[name]["threads"].add( str( p.num_threads( ) ) )
#     for children in p.children( recursive=True ):
#         children_pid = children.pid
#         if psutil.pid_exists( children_pid ):
#             try:
#                 name = "%s####%s###%s" % (str( children_pid ), " ".join( children.cmdline( ) ), str( children.create_time( ) ))
#                 process_info_consumtion.setdefault( name, {} )
#                 process_info_consumtion[name].setdefault( "memory", set( ) )
#                 process_info_consumtion[name]["memory"].add( str( children.memory_info( )[0] ) )
#                 process_info_consumtion[name].setdefault( "threads", set( ) )
#                 process_info_consumtion[name]["threads"].add( str( children.num_threads( ) ) )
#             except:
#                 pass
#
# print(process_info_consumtion)
#
# #CAMITAX
# cmd = shlex.split("kraken --fastq-input --paired --gzip-compressed --threads 5 --db /home/bfosso/share/kraken_db /home/bfosso/simulated_microbiome/upgrade_MetaShot/simulated_metagenome_1.fastq.gz /home/bfosso/simulated_microbiome/upgrade_MetaShot/simulated_metagenome_2.fastq.gz --output symbiome.kraken")
# p = psutil.Popen( cmd )
# process_id = p.pid
# while pid_status( p.pid ) not in ["finished", "zombie"]:
#     name = "%s####%s###%s" % (str( process_id ), " ".join( p.cmdline( ) ), str( p.create_time( ) ))
#     process_info_consumtion.setdefault( name, {} )
#     process_info_consumtion[name].setdefault( "memory", set( ) )
#     process_info_consumtion[name]["memory"].add( str( p.memory_info( )[0] ) )
#     process_info_consumtion[name].setdefault( "threads", set( ) )
#     process_info_consumtion[name]["threads"].add( str( p.num_threads( ) ) )
#     for children in p.children( recursive=True ):
#         children_pid = children.pid
#         if psutil.pid_exists( children_pid ):
#             try:
#                 name = "%s####%s###%s" % (str( children_pid ), " ".join( children.cmdline( ) ), str( children.create_time( ) ))
#                 process_info_consumtion.setdefault( name, {} )
#                 process_info_consumtion[name].setdefault( "memory", set( ) )
#                 process_info_consumtion[name]["memory"].add( str( children.memory_info( )[0] ) )
#                 process_info_consumtion[name].setdefault( "threads", set( ) )
#                 process_info_consumtion[name]["threads"].add( str( children.num_threads( ) ) )
#             except:
#                 pass
#
# # cmd= shlex.split("kraken-translate --db /home/bfosso/share/kraken_db symbiome.kraken --mpa-format > symbiome.labels")
# # p = psutil.Popen( cmd )
# # process_id = p.pid
# # while pid_status( p.pid ) not in ["finished", "zombie"]:
# #     name = "%s####%s###%s" % (str( process_id ), " ".join( p.cmdline( ) ), str( p.create_time( ) ))
# #     process_info_consumtion.setdefault( name, {} )
# #     process_info_consumtion[name].setdefault( "memory", set( ) )
# #     process_info_consumtion[name]["memory"].add( str( p.memory_info( )[0] ) )
# #     process_info_consumtion[name].setdefault( "threads", set( ) )
# #     process_info_consumtion[name]["threads"].add( str( p.num_threads( ) ) )
# #     for children in p.children( recursive=True ):
# #         children_pid = children.pid
# #         if psutil.pid_exists( children_pid ):
# #             try:
# #                 name = "%s####%s###%s" % (str( children_pid ), " ".join( children.cmdline( ) ), str( children.create_time( ) ))
# #                 process_info_consumtion.setdefault( name, {} )
# #                 process_info_consumtion[name].setdefault( "memory", set( ) )
# #                 process_info_consumtion[name]["memory"].add( str( children.memory_info( )[0] ) )
# #                 process_info_consumtion[name].setdefault( "threads", set( ) )
# #                 process_info_consumtion[name]["threads"].add( str( children.num_threads( ) ) )
# #             except:
# #                 pass
#
# #METAPHLAN2
# cmd = shlex.split("/home/bfosso/Codici_programmi/biobakery-metaphlan2-9428f64c98ad/metaphlan2.py simulated_metagenome_1.fastq,simulated_metagenome_2.fastq --input_type fastq --nproc 5 --bowtie2out metagenome.bowtie2.bz2  -t rel_ab -o profiled_metagenome.txt")
# p = psutil.Popen( cmd )
# process_id = p.pid
# while pid_status( p.pid ) not in ["finished", "zombie"]:
#     name = "%s####%s###%s" % (str( process_id ), " ".join( p.cmdline( ) ), str( p.create_time( ) ))
#     process_info_consumtion.setdefault( name, {} )
#     process_info_consumtion[name].setdefault( "memory", set( ) )
#     process_info_consumtion[name]["memory"].add( str( p.memory_info( )[0] ) )
#     process_info_consumtion[name].setdefault( "threads", set( ) )
#     process_info_consumtion[name]["threads"].add( str( p.num_threads( ) ) )
#     for children in p.children( recursive=True ):
#         children_pid = children.pid
#         if psutil.pid_exists( children_pid ):
#             try:
#                 name = "%s####%s###%s" % (str( children_pid ), " ".join( children.cmdline( ) ), str( children.create_time( ) ))
#                 process_info_consumtion.setdefault( name, {} )
#                 process_info_consumtion[name].setdefault( "memory", set( ) )
#                 process_info_consumtion[name]["memory"].add( str( children.memory_info( )[0] ) )
#                 process_info_consumtion[name].setdefault( "threads", set( ) )
#                 process_info_consumtion[name]["threads"].add( str( children.num_threads( ) ) )
#             except:
#                 pass
#
# script_time = time.strftime( "%d/%m/%Y %H:%M:%S", time.localtime( time.time( ) ) )
# print("END %s" % script_time)
#
# tmp = open( "technical_measurement.tsv", "w" )
# for process in process_info_consumtion.keys( ):
#     tmp.write( "%s\t%s\t%s\n" % (process, ";".join( list( process_info_consumtion[process]["memory"] ) ), ";".join( list( process_info_consumtion[process]["threads"] ) )) )
# tmp.close( )
