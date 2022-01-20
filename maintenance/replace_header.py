#
# Copyright 2022 Johannes Laurin Hoermann
#
# ### MIT license
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#


import sys
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

file_lines = open(sys.argv[1], 'r').readlines()
header_lines = sys.stdin.readlines()

# keep jupytext header at top if present
while file_lines[0].startswith('#') and not file_lines[0].startswith('# ---'):
    file_lines = file_lines[1:]

top_comment_lines = []
# This block aims at Jupytext files and expects two '# ---' markers
if file_lines[0].startswith('# ---'):
    logger.debug("Found '# ---#' top comment start in file.")
    # mimic do ... while
    while file_lines[0].startswith('#'):
        top_comment_lines.append(file_lines[0])
        file_lines = file_lines[1:]
        if file_lines[0].startswith('# ---'):
            logger.debug("Found '# ---#' top comment end in file.")
            logger.debug("Kept '{} top comment.".format(top_comment_lines))
            break
    top_comment_lines.append(file_lines[0])

while file_lines[0].startswith('#'):
    file_lines = file_lines[1:]

file_lines.insert(0, '#\n')
for header_line in header_lines[::-1]:
    file_lines.insert(0, '# {}'.format(header_line).strip() + '\n')
file_lines.insert(0, '#\n')

for top_comment_line in top_comment_lines[::-1]:
    file_lines.insert(0, (top_comment_line).strip() + '\n')

open(sys.argv[1], 'w').writelines(file_lines)
