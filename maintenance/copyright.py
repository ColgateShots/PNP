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


import os
import sys
import logging
from collections import defaultdict
from datetime import datetime
from subprocess import Popen, PIPE

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

root = os.path.dirname(sys.argv[0])


def read_authors(fn):
    return {email.strip('<>'): name for name, email in
            [line.rsplit(maxsplit=1) for line in open(fn, 'r')]}


def merge_committers(lhs_committers, rhs_committers):
    """Union between lhs_ and rhs_committers, returns modified lhs_committers"""
    for committer, years in rhs_committers.items():
        if committer in lhs_committers:
            lhs_committers[committer] |= years
        else:
            lhs_committers[committer] = years
    return lhs_committers


def parse_git_log(log, authors):
    committers = defaultdict(set)
    author = None
    date = None
    for line in log.decode('utf-8').split('\n'):
        logger.debug("Parse '{}'".format(line))
        if line.startswith('commit'):
            logger.debug("Parse commit")
            if date is not None and author is not None:
                committers[author].add(date.year)
                logger.debug("Added year of date '{}' to author '{}'.".format(author, date))
        elif line.startswith('Author:'):
            email = line.rsplit('<', maxsplit=1)[1][:-1]
            logger.debug("Parsed email '{}'.".format(email))
        elif line.startswith('Date:'):
            date = datetime.strptime(line[5:].rsplit(maxsplit=1)[0].strip(),
                                     '%a %b %d %H:%M:%S %Y')
            logger.debug("Parsed date '{}'.".format(date))

            try:
                logger.debug("Try to get author name from eMail '{}'.".format(email))
                author = authors[email]
            except KeyError:
                logger.debug("No entry for author, use email '{}'.".format(email))
                author = email
        # elif 'copyright' in line.lower() or 'license' in line.lower():
        #    date = None
        # why should this exception be necessary? It results in erroneously ignored commits
    logger.debug("Done parsing bulk. Last parset author: {}'.".format(author))
    if date is not None:
        committers[author].add(date.year)
        logger.debug("Added year of date '{}' to final author entry.".format(date))
    return committers


def pretty_years(years):
    def add_to_year_string(s, pprev_year, prev_year):
        if pprev_year == prev_year:
            # It is a single year
            if s is None:
                return f'{prev_year}'
            else:
                return f'{s}, {prev_year}'
        else:
            # It is a range
            if s is None:
                return f'{pprev_year}-{prev_year}'
            else:
                return f'{s}, {pprev_year}-{prev_year}'

    years = sorted(years)
    prev_year = pprev_year = years[0]
    s = None
    for year in years[1:]:
        if year - prev_year > 1:
            s = add_to_year_string(s, pprev_year, prev_year)
            pprev_year = year
        prev_year = year
    return add_to_year_string(s, pprev_year, prev_year)


authors = read_authors('{}/../AUTHORS'.format(root))
logger.info("Read {} from AUTHORS file".format(authors))

# if repository contains jupytext-synced ipynb files, also parse their authors
py_filename = sys.argv[1]
basename, extension = os.path.splitext(py_filename)
ipynb_filename = basename  + '.ipynb'

process = Popen(['git', 'log', '--follow', py_filename], stdout=PIPE,
                stderr=PIPE)
stdout, stderr = process.communicate()
logger.info("Retrieved '{}' from 'git log --follow {}'.".format(stdout, py_filename))
committers = parse_git_log(stdout, authors)
logger.info("Parsed {} from {}.".format(committers, py_filename))

if os.path.isfile(ipynb_filename):
    logger.info("Found according ipynb notebookt {}.".format(ipynb_filename))
    process = Popen(['git', 'log', '--follow', ipynb_filename], stdout=PIPE,
                    stderr=PIPE)
    stdout, stderr = process.communicate()
    logger.info("Retrieved '{}' from 'git log --follow {}'.".format(stdout, ipynb_filename))
    ipynb_committers = parse_git_log(stdout, authors)
    logger.info("Parsed {} from {}.".format(ipynb_committers, ipynb_filename))
    committers = merge_committers(committers, ipynb_committers)
    logger.info("Merged committers into {}.".format(committers))

prefix = 'Copyright'
for name, years in committers.items():
    print('{} {} {}'.format(prefix, pretty_years(years), name))
    prefix = ' ' * len(prefix)
print()
