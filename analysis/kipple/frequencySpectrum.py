#! /usr/bin/env python2.7

from HelPopGen import dataRecords, phaseRecord
from collections import Counter

def getFreqs(record):
    ls = len(record[0][1]) - len(record[0][1]) % 3
    sites = [[x[1][i:i+3] for x in record] for i in range(0, ls, 3)]
    usable_sites = [s for s in sites if
            not any(["N" in c for c in s]) and
            not any([c in ["TAG","TGA","TAA"] for c in s])
            and len(set(s)) <= 2] # include only biallelic and constant sites
    freqs_list = [min(Counter(s).values()) for s in usable_sites]
    return freqs_list

# frequency spectrum for aglaope and amaryllis together
phased_records = (phaseRecord(record) for record in dataRecords())
agam = ([x for x in record if x[0].startswith("a")] for record in phased_records)
freqs = list()
for record in agam:
    freqs.extend(getFreqs(record))
frequency_spectrum = Counter(freqs)
frequency_spectrum[0] = frequency_spectrum[max(frequency_spectrum.keys())]
del(frequency_spectrum[max(frequency_spectrum.keys())])
print "amaryllis and aglaope:  ", frequency_spectrum
# amaryllis and aglaope:   Counter({0: 4852130, 1: 254392, 2: 81829, 3: 46117, 
#         4: 32798, 5: 27416, 6: 23848, 7: 22502, 8: 12977})

# frequency spectrum for amaryllis only
phased_records = (phaseRecord(record) for record in dataRecords())
aglaope = ([x for x in record if x[0].startswith("amaryllis")] for record in phased_records)
freqs = list()
for record in aglaope:
    freqs.extend(getFreqs(record))
frequency_spectrum = Counter(freqs)
frequency_spectrum[0] = frequency_spectrum[max(frequency_spectrum.keys())]
del(frequency_spectrum[max(frequency_spectrum.keys())])
print "amaryllis:  ", frequency_spectrum
# amaryllis:   Counter({0: 5084918, 1: 208634, 2: 78070, 3: 53729, 4: 26841})

# frequency spectrum for aglaope only
phased_records = (phaseRecord(record) for record in dataRecords())
aglaope = ([x for x in record if x[0].startswith("aglaope")] for record in phased_records)
freqs = list()
for record in aglaope:
    freqs.extend(getFreqs(record))
frequency_spectrum = Counter(freqs)
frequency_spectrum[0] = frequency_spectrum[max(frequency_spectrum.keys())]
del(frequency_spectrum[max(frequency_spectrum.keys())])
print "aglaope:  ", frequency_spectrum
# aglaope:   Counter({0: 5074911, 1: 214229, 2: 77164, 3: 53341, 4: 27363})

# frequency spectrum for erato only
phased_records = (phaseRecord(record) for record in dataRecords())
aglaope = ([x for x in record if x[0].startswith("erato")] for record in phased_records)
freqs = list()
for record in aglaope:
    freqs.extend(getFreqs(record))
frequency_spectrum = Counter(freqs)
frequency_spectrum[0] = frequency_spectrum[max(frequency_spectrum.keys())]
del(frequency_spectrum[max(frequency_spectrum.keys())])
print "erato:  ", frequency_spectrum
# erato:   Counter({0: 3920087, 1: 231479, 2: 92935, 3: 61291, 4: 32726})

# count sites that segregate between aglaope and amaryllis
phased_records = (phaseRecord(record) for record in dataRecords())
agam = ([x for x in record if x[0].startswith("a")] for record in phased_records)
count = 0
for record in agam:
    ls = len(record[0][1]) - len(record[0][1]) % 3
    sites = [[x[1][i:i+3] for x in record] for i in range(0, ls, 3)]
    for s in sites:
        if not any(["N" in c for c in s]):
            if not any([c in ["TGA","TAG","TAA"] for c in s]):
                if len(set(s)) == 2:
                    if len(set(zip(s, [x[0][:2] for x in record]))) == 2:
                        count += 1
print "aglaope/amaryllis segregating sites:  ", count
# The result is 46 which is actually more than we would expect by chance, so good.

# Count sites that are at 50% frequency and are heterozygous in 8 individuals
phased_records = (phaseRecord(record) for record in dataRecords())
agam = ([x for x in record if x[0].startswith("a")] for record in phased_records)
count = 0
for record in agam:
    ls = len(record[0][1]) - len(record[0][1]) % 3
    sites = [[x[1][i:i+3] for x in record] for i in range(0, ls, 3)]
    for s in sites:
        if not any(["N" in c for c in s]):
            if not any([c in ["TGA","TAG","TAA"] for c in s]):
                if len(set(s)) == 2:
                    if all([s[i] != s[i+1] for i in range(0, len(s), 2)]):    
                        count += 1
print "aglaope/amaryllis segregating sites:  ", count
# The result is 2628

