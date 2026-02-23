#  MIT License
#
#  Copyright (c) 2021-2026. Aleksandr Serdiukov, Anton Zamyatin, Aleksandr Sinitsyn, Vitalii Dravgelis and Computer Technologies Laboratory ITMO University team.
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.

from typing import List
from setuptools import find_packages, setup


requirements: List[str] = []
with open("requirements.txt", mode="rt", encoding="utf-8") as f:
    requirements = f.readlines()

setup(
    name='hict_server',
    version='0.1.3rc1',
    packages=list(set(['hict_server', 'hict_server.api_controller',
              'hict_server.api_controller.dto']).union(find_packages())),
    url='https://genome.ifmo.ru',
    license='MIT',
    author='Alexander Serdiukov, Anton Zamyatin and CT Lab ITMO University team',
    author_email='',
    description='Development version of API and tiling server for HiCT interactive Hi-C scaffolding tool.',
    install_requires=list(set([
        'hict>=0.1.3rc1,<1.0',
        'hict_utils>=0.1.3rc1,<1.0',
    ]).union(requirements)),
    entry_points={
        'console_scripts': ['hict_server=hict_server.api_controller.dev_demo_server:main'],
    }
)
