from setuptools import setup, Extension
import os

with open(os.path.join(os.path.dirname(__file__), "aireal", "version.py")) as f_in:
    exec(f_in.read())

package =  {"name": "pipeline",
            "version": __version__,
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
                                                  "call_variants=pipeline.call_variants:main",
                                                  "filter_vcf=pipeline.filter_vcf:main",
                                                  "postprocess_mutect2_vcf=pipeline.postprocess_mutect2_vcf:main",
                                                  "postprocess_varscan_vcf=pipeline.postprocess_varscan_vcf:main",
                                                  "multiplexing=pipeline.multiplexing:main",
                                                  "mount_instance_storage=pipeline.aws:mount_instance_storage",
                                                  "bscopy=pipeline.support.bscopy:main",
                                                  "profile=pipeline.support.profile:main",
                                                  "s3_wrap=pipeline.support.s3_wrap:main",
                                                  "sqs_dequeue=pipeline.support.sqs_dequeue:main",
                                                  "sqs_enqueue=pipeline.support.sqs_enqueue:main",
                                                  "chr_rename_fasta=pipeline.support.chr_rename_fasta:main",
                                                  "covermi_stats=pipeline.covermi_stats:main",
                                                  "vcf_stats=pipeline.vcf_stats:main",
                                                  "ontarget=pipeline.ontarget:main",
                                                  "annotate_panel=pipeline.annotate_panel:main",
                                                  "size=pipeline.size:main",
                                                  "breakpoint=pipeline.breakpoint:main"]},
            }

setup(**package)









