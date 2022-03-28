

def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b:
            break
        yield b


def main(filename):
    print('Counting number of lines in file: %s' % filename)
    # count number of lines in sam file
    with open(filename, "r", encoding="utf-8", errors='ignore') as f:
        linecount = sum(bl.count("\n") for bl in blocks(f))

    print(linecount)

    return linecount