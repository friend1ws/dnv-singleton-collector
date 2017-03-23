#!/usr/bin/env python

from distutils.core import setup

setup(name='dnv_singleton_collector',
      version='0.1.0',
      description='Script for collecting DNV singleton from bam files',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/friend1ws/dnv_singleton_collector.git',
      package_dir = {'': 'lib'},
      packages=['dnv_singleton_collector'],
      scripts=['dnv_singleton_collector'],
      license='GPL-3'
     )

