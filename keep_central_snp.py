#!/usr/bin/env python3

# Copyright 2019 Duarte Teomoteo Balata <duarte.balata@gmail.com>
# keep_central_snps.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# keep_central_snps is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with keep_central_snps. If not, see <http://www.gnu.org/licenses/>.

import argparse
import sys

def get_args(args):
    parser = argparse.ArgumentParser(prog='python3')
    
    parser.add_argument("input", metavar="input.vfc", type=str,
                        help="VCF file with all the SNPs.")
    
    parser.add_argument("-o","--output", metavar="output.vcf", dest="output", type=str,
                        help="output file with the SNPs closest to the centre of each locus")
    
    parser.add_argument("-l","--len", metavar="int", dest="length", type=str,
                        help="length of each locus")
    
    arguments = parser.parse_args(args)
    return arguments

def write_vcf_headers(vcf_path,output_path):
    vcf_file= open(vcf_path,"r")
    out_file=open(output_path,"w")
    
    for line in vcf_file:
        if line.startswith("#"):
            out_file.write(line)
        else:
            break
    
    vcf_file.close()
    out_file.close()
    

def write_vcf_body(vcf_path,output_path,locus_length):
    
    vcf_file= open(vcf_path,"r")
    out_file=open(output_path,"a")
    
    buff=[]
    stop_position=0
    last_dist=0
    dist=0
    last_name=""
    
    for line in vcf_file:    
        try:  
            name=line.split("\t")[0]
            dist = int(line.split("\t")[1])
            if dist >= stop_position or dist < last_dist or last_name != name:
                
                if len(buff)>0:
                    out_file.write(buff[len(buff)//2])
                    
                buff=[]
                stop_position = dist + int(locus_length)
                buff.append(line)
            else:
                buff.append(line)
            last_dist=dist
            last_name=name
        except:
            pass
    
    out_file.write(buff[len(buff)//2])
        
    vcf_file.close()        
    out_file.close()

if __name__ == "__main__":
    args= get_args(sys.argv[1:])
    write_vcf_headers(args.input,args.output)
    write_vcf_body(args.input,args.output,args.length)     
