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
                                                  "mount_instance_storage=pipeline.aws:mount_instance_storage",
                                                  "udini=pipeline.udini:main",
                                                  "elduderino=pipeline.elduderino:main",
                                                  "vcf_pvalue_2_phred=pipeline.vcf_pvalue_2_phred:main",
                                                  "covermi_stats=pipeline.covermi_stats:main",
                                                  "vcf_stats=pipeline.vcf_stats:main",
                                                  "ontarget=pipeline.ontarget:main",
                                                  "annotate_panel=pipeline.annotate_panel:main",
                                                  "size=pipeline.size:main",
                                                  "trim_sam=pipeline.trim_sam:main",
                                                  "breakpoint=pipeline.breakpoint:main",
                                                  "bscopy=pipeline.bscopy:main",
                                                  "enqueue=pipeline.enqueue:main",
                                                  "wrap_task=pipeline.wrap_task:main",
                                                  "runner=pipeline.runner:main"] },
            }

setup(**package)









