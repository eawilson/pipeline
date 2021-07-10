from setuptools import setup, Extension

package =  {"name": "pipeline",
            "version": "0.1",
            "description": "NGS sequencing bioinformatics pipeline",
            "url": "",
            "author": "Ed Wilson",
            "author_email": "edwardadrianwilson@yahoo.co.uk",
            "license": "MIT",
            "packages": ["pipeline"],
            "install_requires": ["covermi"],
            "include_package_data": True,
            "zip_safe": True,
            "entry_points": { "console_scripts": ["cfpipeline=pipeline.cfpipeline:main",
                                                  "cfpipeline2=pipeline.cfpipeline2:main",
                                                  "call_variants=pipeline.call_variants:main",
                                                  "multiplexing=pipeline.multiplexing:main",
                                                  "mount_instance_storage=pipeline.aws:mount_instance_storage",
                                                  "bscopy=pipeline.support.bscopy:main",
                                                  "enqueue=pipeline.support.enqueue:main",
                                                  "runner=pipeline.support.runner:main",
                                                  "profile=pipeline.support.profile:main",
                                                  "s3_wrap=pipeline.support.s3_wrap:main",
                                                  "sqs_dequeue=pipeline.support.sqs_dequeue:main",
                                                  "sqs_enqueue=pipeline.support.sqs_enqueue:main",
                                                  "chr_rename_fasta=pipeline.support.chr_rename_fasta:main",
                                                  "udini=pipeline.udini:main",
                                                  "elduderino_py=pipeline.elduderino_py:main",
                                                  "vcf_pvalue_2_phred=pipeline.vcf_pvalue_2_phred:main",
                                                  "covermi_stats=pipeline.covermi_stats:main",
                                                  "vcf_stats=pipeline.vcf_stats:main",
                                                  "ontarget=pipeline.ontarget:main",
                                                  "annotate_panel=pipeline.annotate_panel:main",
                                                  "size=pipeline.size:main",
                                                  "trim_sam=pipeline.trim_sam:main",
                                                  "breakpoint=pipeline.breakpoint:main"]},
            }

setup(**package)









