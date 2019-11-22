from setuptools import setup, Extension

package =  {"name": "pipeline",
            "version": "0.1",
            "description": "NGS sequencing bioinformatics pipeline",
            "url": "",
            "author": "Ed Wilson",
            "author_email": "edwardadrianwilson@yahoo.co.uk",
            "license": "MIT",
            "packages": ["pipeline"],
            "ext_modules": [Extension("pipeline.cpipeline", ["pipeline/cpipeline/pipeline.c"])],
            "install_requires": [],
            "include_package_data": True,
            "zip_safe": True,
            "entry_points": { "console_scripts": ["cfpipeline=pipeline.cfpipeline:main"] },
            }

setup(**package)









