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
            "entry_points": { "console_scripts": ["cfpipeline=pipeline.pipelines.cfpipeline:main",
                                                  "filter_sam=pipeline.filter_sam:main",
                                                  "terminate=pipeline.runner:terminate",
                                                  "udini=pipeline.udini:main",
                                                  "elduderino=pipeline.elduderino:main",
                                                  "vcf_pvalue_2_phred=pipeline.vcf_pvalue_2_phred:main",
                                                  "runner=pipeline.runner:main"] },
            }

setup(**package)









