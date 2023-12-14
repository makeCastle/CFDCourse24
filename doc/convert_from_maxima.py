import sys
import re

re_power = r'((\w|\d)+)\^(\d)'
re_power_query = re.compile(re_power, re.M)

def convert_string(s):
    # remove square parenthesis
    s = s.replace('[', '')
    s = s.replace(']', '')
    # powers
    for it in reversed([x for x in re_power_query.finditer(s)]):
        arg = s[it.start(1):it.end(1)]
        pwr = int(s[it.start(3):it.end(3)])
        changer = arg
        for i in range(2, pwr+1):
            changer += "*" + arg
        s = s[:it.start(0)] + changer + s[it.end(0):]
    return s
    
txtin = []

for s in sys.stdin:
    s = s.strip()
    if not s:
        break
    txtin.append(s)


print('\n\n\n')

for s in txtin:
    print(convert_string(s))
