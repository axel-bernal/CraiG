import sys, re

class CraigUndefined(Exception):
    pass

def readRecords(inputFile,
                recordSep="\\n",
                readSize=8192):
    """Like the normal file iter but you can set what string indicates newline.

    The newline string can be arbitrarily long; it need not be restricted to a
    single character. You can also set the read size and control whether or not
    the newline string is left on the end of the iterated lines.  Setting
    newline to '\0' is particularly good for use with an input file created with
    something like "os.popen('find -print0')".
    """

    partialLine = ''
    mrs = readSize
    while True:
        charsJustRead = ''
        chunk = inputFile.read(mrs)
        # maybe record is larger than default size..

        while chunk:
            charsJustRead += chunk
            if charsJustRead.find(recordSep) >= 0: break
            mrs = 2*mrs
            chunk = inputFile.read(mrs)

        if not charsJustRead:
            break
        partialLine += charsJustRead
        lines = partialLine.split(recordSep)
        partialLine = lines.pop()
        for line in lines:
            yield line
    if partialLine: yield partialLine


recParsei = re.compile(r"^\>?(\S+)(\s*[\S,\s,\n]*)$")
recParsec = re.compile(r"^\>?(\S+)(\s+\d?\<?)([^\;,\],\[,\<,\s]+)(\_[^\_,\],\[,\;]+\_[^\_,\],\[,\;]+.*)$")
recParsel = re.compile(r"^\>?(\S+)(\s+\d?\<?)([^\s,\>]+)(.*)$")

def processRecord(line, idType):
    if idType == 'i':
        match = recParsei.match(line)
        if match:
            return (match.group(1).lstrip('>'), match.group(2).rstrip('\n'))
    elif idType == 'c':
        match = recParsec.match(line)
        if match:
            return (match.group(3), line.rstrip('\n').lstrip('>'))
    elif idType == 'l':
        match = recParsel.match(line)
        if match:
            return (match.group(3), line.rstrip('\n').lstrip('>'))
    else: return (None, None)

def printRecord(key, value, recordSep, idType):
    index_nl = recordSep.find('\n')
    if idType == 'i':
        print recordSep[index_nl + 3:] + key + \
        value + (recordSep[:index_nl] if index_nl > 0 else "")
    elif idType == 'c' or idType == 'l':
        print recordSep[index_nl + 3:] + value + \
        (recordSep[:index_nl] if index_nl > 0 else "")

exonParse = re.compile(r"^\s*(\S+)_(\d+)_(\d+)\s*")


def readLocExons(location, exons) :
    exons = []
    res = location.split(";")
    for i in range(len(res)):
        match = exonParse.match(res[i])
        if match:
            exons.append([match.group(1),match.group(2),match.group(3)])
        else:
            raise CraigUndefined("Wrong exon format in "+location)


wlocParse = re.compile(r"^([^\s]+)\s+([\d,\.,\,\-;]+)$")
fplocParse = re.compile(r"([^\[]+)\[(\S+)\s*$")
tplocParse = re.compile(r"([^\]]+)\](\S+)\s*$")
locParse1 = re.compile(r"(\d+)?\<(\S*)\s*$")
locParse2 = re.compile(r"^(\S*)\>/")
locParse3 = re.compile(r"^([^\s]+)$")

def readLocation(location, fiveUtrExons, exons, threeUtrExons):
    fiveUtrExons = []
    exons = []
    threeUtrExons = []
    noStart = 0
    phase = 0
    noStop = 0
    weight = []
    match = wlocParse.match(location)

    if match: # weighted locations
	location = match.group(1)
	weight = match.group(2).split(";")

    match = fplocParse.match(location)
    if match:
        location = match.group(2)
        readLocExons(match.group(1), fiveUtrExons)

    match = tplocParse.match(location)
    if match:
        location = match.group(1)
        readLocExons(match.group(2), threeUtrExons)

    match = locParse1.match(location)
    if match :
        phase = match.group(1)
        location = match.group(2)
	noStart = 1 # partial first exon
#    print location, match
    match = locParse2.match(location)
    if match :
        location = match.group(1)
        noStop = 1 # partial last exon

    match = locParse3.match(location)
    if match :
        readLocExons(location, exons)
        for i in range(len(weight)):
	    exons[i].append(weight[i])
    else :
        raise CraigUndefined("Wrong exon format in "+location)
