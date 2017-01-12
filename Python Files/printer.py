
## returns the input with a timestamp

import time

def printer(input):
    now = time.strftime("%c")
    output = input + ' - ' + now
    return output
