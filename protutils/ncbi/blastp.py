"""
Taken from ProDy
(http://www.csb.pitt.edu/prody/_modules/prody/proteins/blastpdb.html)
"""
import re
import time
import urllib2
import xml.etree.cElementTree as etree
from urllib import urlencode


def blast_pdb(sequence, nhits=250, expect=1e-10, timeout=60, pause=1):
    query = {
        'DATABASE': 'pdb',
        'ENTREZ_QUERY': '(none)',
        'PROGRAM': 'blastp',
        'EXPECT': expect,
        'HITLIST_SIZE': nhits,
        'CMD': 'Put',
        'QUERY': sequence
    }

    url = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi'
    data = urlencode(query)

    request = urllib2.Request(
        url, data=data, headers={'User-agent': 'protutils'}
    )
    response = urllib2.urlopen(request)

    html = response.read()
    m = re.search('RID =\s?(.*?)\n', html)
    if m:
        rid = m.group(1)
    else:
        raise Exception('Could not parse response.')

    query = {
        'ALIGNMENTS': 500,
        'DESCRIPTIONS': 500,
        'FORMAT_TYPE': 'XML',
        'RID': rid,
        'CMD': 'Get'
    }
    data = urlencode(query)

    slept = 0
    while slept < timeout:
        request = urllib2.Request(
            url, data=data, headers={'User-agent': 'protutils'}
        )
        response = urllib2.urlopen(request)
        results = response.read()
        m = re.search('Status=(.*?)\n', results)
        if not m:
            break
        elif m.group(1).strip().upper() == 'READY':
            break
        else:
            time.sleep(pause)
            slept += pause
    with open('blastp.xml', 'w') as f:
        f.write(results)
    return etree.XML(results)


def xml_dict(root, tag_prefix):
    d = {}
    regex = re.compile(r'{0}(.*)'.format(tag_prefix))

    for element in root:
        tag = element.tag
        m = regex.search(tag)
        if m:
            key = m.group(1)
            if len(element) == 0:
                d[key] = element.text
            else:
                d[key] = element
    return d


class BLASTPDBRecord(object):

    def __init__(self, sequence, nhits=250, expect=1e-10, timeout=60, pause=1):
        self.qseq = sequence
        root = blast_pdb(sequence, nhits, expect, timeout, pause)
        root = xml_dict(root, 'BlastOutput_')
        self.query_id = root['query-ID']
        if not len(sequence) == int(root['query-len']):
            raise ValueError('Sequence length does not match query length')
        self.param = xml_dict(root['param'][0], 'Parameters_')

        hits = []
        for elem in root['iterations']:
            for child in xml_dict(elem, 'Iteration_')['hits']:
                hit = xml_dict(child, 'Hit_')
                data = xml_dict(hit['hsps'][0], 'Hsp_')
                for key in ['align-len', 'gaps', 'hit-frame', 'hit-from',
                            'hit-to', 'identity', 'positive', 'query-frame',
                            'query-from', 'query-to']:
                    data[key] = int(data[key])
                for key in ['evalue', 'bit-score', 'score']:
                    data[key] = float(data[key])
                p_identity = (data['identity'] /
                              float(data['query-to'] - data['query-from'] + 1)
                              * 100)
                p_overlap = ((data['align-len'] - data['gaps']) /
                             float(len(sequence)) * 100)
                data['percent_identity'] = p_identity
                data['percent_overlap'] = p_overlap
                __, gi, __, pdb, chain = hit['id'].split('|')
                data['gi'] = gi
                data['pdb'] = pdb
                data['chain'] = chain
                data['def'] = hit['def']
                hits.append(data)
        hits.sort(key=lambda x: x['percent_identity'], reverse=True)
        self.hits = hits

    def get_hits(self, percent_identity=90.0, percent_overlap=70.0):
        hits = {}
        for hit in self.hits:
            if hit['percent_identity'] < percent_identity:
                break
            if hit['percent_overlap'] < percent_overlap:
                continue
            key = '{pdb}_{chain}'.format(**hit)
            hits[key] = hit
        return hits

    def get_best(self):
        return self.hits[0]

    def ranking(self):
        return {
            '{pdb}_{chain}'.format(**hit): hit[
                'percent_identity'
            ] for hit in self.hits
        }
