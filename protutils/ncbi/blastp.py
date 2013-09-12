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
