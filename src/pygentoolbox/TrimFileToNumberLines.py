def main(f, limit, gz=False):
    # f = str, limit = int
    # f is full path to file # limit is line number limit 1st line is 1 NOT 0 # still outputs the linecount <= limit
    if limit < 1000:
        outf = '.'.join(f.split('.')[:-1] + ['lines%d' % limit, f.split('.')[-1]])
    elif limit < 1000000:
        outf = '.'.join(f.split('.')[:-1] + [f'lines{round(limit/1000)}k', f.split('.')[-1]])
    else:
        outf = '.'.join(f.split('.')[:-1] + [f'lines{round(limit/1000000)}M', f.split('.')[-1]])

    if gz is False:
        with open(outf, 'w') as OUT:
            OUT.write('')
        checkpoint = int(limit / 4)
        with open(f, 'r') as FILE:
            with open(outf, 'a') as OUT:
                for count, line in enumerate(FILE, start=1):
                    if count <= limit:
                        OUT.write(line)
                        if count % checkpoint == 0:
                            print('On line: %d' % count)
                            if count == limit:
                                break
    elif gz is True:
        import gzip
        with gzip.open(outf, 'wb') as OUT:
            OUT.write(b'')
        checkpoint = int(limit / 4)
        with gzip.open(f, 'rb') as FILE:
            with gzip.open(outf, 'ab') as OUT:
                for count, line in enumerate(FILE, start=1):
                    if count <= limit:
                        OUT.write(line)
                        if count % checkpoint == 0:
                            print('On line: %d' % count)
                            if count == limit:
                                break

    print(outf)
    print('Finished')
