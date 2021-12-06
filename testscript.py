import argparse

def main(args):
    with open(args.output, 'w') as fo, open(args.input, 'r') as fi:
        modules = []
        for line in fi:
            u, v, _ = line.split('\t')
            added = False
            for module in modules:
                if u in module or v in module:
                    module |= {u, v}
                    added = True
            if not added:
                modules.append({u, v})
        for i, module in enumerate(modules):
            fo.write('\n'.join(['{}\t{}'.format(v, i) for v in module]) + '\n')


if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type=str, required=True,
                        help="input graph file name")
    parser.add_argument('-o', '--output', type=str, required=True, 
                        help="output module file name")

    main(parser.parse_args())
