
# def oshdufoshdfo_web_logo(fastafile, title):
#     import weblogo
#     import io
#     import PIL.Image
#     import os
#
#     # args = FILE
#
#     infile = open(fastafile)
#     seqs = weblogo.read_seq_data(infile)
#     logodata = weblogo.LogoData.from_seqs(seqs)
#     logooptions = weblogo.LogoOptions()
#     logooptions.title = title
#     logoformat = weblogo.LogoFormat(logodata, logooptions)
#     eps = weblogo.eps_formatter(logodata, logoformat)
#     png = weblogo.png_print_formatter(logodata, logoformat)
#
#     image = PIL.Image.open(io.BytesIO(png))
#     path, f = os.path.split(fastafile)
#     outfile = os.path.join(path, f'{title}.png')
#     image.save('IMAGE.PNG')


def web_logo(fastafile, title, molecule, outformat, extra):
    import os
    path, f = os.path.split(fastafile)
    outfile = os.path.join(path, f'{title}.{outformat}')
    # weblogo -f in.fa -t title -F pdf -A dna > outfile
    args = f'-f {fastafile} -t {title} -A {molecule} -F {outformat} {extra}'
    cmd = f'weblogo {args} > {outfile}'
    print(cmd)
    os.system(cmd)

    return outfile


def main(fastafile, title, molecule, outformat, extra):
    outfile = web_logo(fastafile, title, molecule, outformat, extra)
    print(f'Weblogo made, output to file: {outfile}')
