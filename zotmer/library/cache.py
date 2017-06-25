import os
import os.path
import sys

import requests

from pykmer.file import openFile

# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000962.3&rettype=fasta&retmode=text

class entrez(object):
    def __init__(self):
        self.base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        self.root = '.zot/cache'
        self.compression = ''

    def __getitem__(self, acc):
        pth = self.makePath(acc, self.compression)

        if not os.path.isfile(pth):
            (d, f) = self.makePathComponents(acc, self.compression)
            if not os.path.isdir(d):
                os.makedirs(d)
            self.download(acc, pth)

        return openFile(pth)

    def download(self, acc, path):
        print >> sys.stderr, 'downloading %s -> %s' % (acc, path)

        qry = {}
        qry['db'] = 'nuccore'
        qry['rettype'] = 'fasta'
        qry['retmode'] = 'text'
        qry['id'] = acc

        # Max size to slurp
        MB64 = 64 * 1024 * 0124

        with requests.get(self.base, params=qry, stream=True) as r:

            # deal with HTTP errors
            r.raise_for_status()

            # write "small" files in one go.
            if 'Content-Length' in r.headers and int(r.headers['Content-Length']) <= MB64:
                with open(path, 'w') as f:
                    f.write(r.content)
                return

            # grab "big" files a piece at a time.
            with open(path, 'w') as f:
                for chk in r.iter_content(chunk_size=None):
                    f.write(chk)

    def makePath(self, acc, compression=''):
        h = hash(acc) & 0xFF
        return '%s/%02x/%s.fasta%s' % (self.root, h, acc, compression)

    def makePathComponents(self, acc, compression=''):
        h = hash(acc) & 0xFF
        return ('%s/%02x' % (self.root, h), '%s.fasta%s' % (acc, compression))

