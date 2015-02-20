def print_val(header):
    line = "Description\tAccession\t"
    header[0]= header[0].split("\t")[-1]
    print header[0]
    print header[-1].rstrip()
    print len(header)
    for i in xrange(len(header)):
        if i < len(header) - 1:
            line += "%s\t" %header[i]
        else:
            line += header[i]
    return line.rstrip()

filename = "/home/ubuntu/ccle/data/CCLE_Expression_2012-09-29.res"
out_file = "array_out.txt"
fp = open(filename, "r")
f_out = open(out_file, "w")
count = 0
for line in fp:
    if count == 0:
        header = line.split("\t\t")
        header = print_val(header)
        print len(header.split("\t"))
        f_out.write(header)
        f_out.write("\n")
    else:
        if count > 0:
            val = ""
            line = line.split("\t")
            if count == 1:
                print line[0], line[1], line[-1]
            #print len(line)
            #print "Length of line: %s" %len(line)
            for i in xrange(len(line) - 1):
                #print i
                if i%2 == 0 or i == 1:

                    if i < (len(line) - 2):
                        val += "%s\t" %line[i]
                    else:
                        val += line[i]
            print len(val.split("\t"))
            f_out.write(val.rstrip())
            f_out.write("\n")

    count += 1
