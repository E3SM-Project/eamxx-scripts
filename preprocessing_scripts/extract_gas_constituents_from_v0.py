#!/usr/bin/env python3
from subprocess import call


def handle_error(return_code):
    if return_code != 0: raise RuntimeError(f'Command failed.')


def main(inputfile, outputfile):

    gas_constituents = ('co', 'ch4', 'co2', 'n2', 'n2o', 'o2', 'o3')
    cmd = f'ncks -5 -O --no-abc -v ' + ','.join(gas_constituents) + f' {inputfile} {outputfile}'
    print(cmd)
    handle_error(call(cmd.split(' ')))


if __name__ == '__main__':
    import plac; plac.call(main)
