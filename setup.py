#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name='minibarcoder',
	description="a pipeline for obtaining barcodes from MinION raw reads",
	version='2.0',
	author='Amrita Srivathsan',
	author_email='asrivathsan@gmail.com',
	packages=find_packages(),
	#packages=['minibarcoder',],
	scripts=['assess_corrbarcodes_wref.py','assess_uncorrbarcodes_wref.py','repool_by_plate.py','consolidate.py','mb_parallel_consensus.py','mb_parallel_demultiplex.py','racon_consensus.py','filter_by_Ns.py','aacorrection.py'],
	package_data={'minibarcoder': ['demfile_2','gplv3.txt','test.fasta']},
	install_requires=['numpy' ,'xlwt'],
	#entry_points={'console_scripts': ['django-admin = django.core.management:execute_from_command_line',]},
	include_package_data=True,
	#data_files=['demfile_2','gplv3.txt','test.fasta'],
	py_modules=['correction_funcs'],
	#long_description=open('README.txt').read(),
	license='GPL v3',
	url='https://github.com/asrivathsan/minibarcoder',
	)
