#################################################################################
# Copyright (c) 2016-, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################
# author: Bo Han (bhan@pacb.com)
# last modified: 2016-07-28 07:48

#########
# Basic #
#########
declare -r MODULE_NAME=Split
declare -r MODULE_VERSION=0.0.2

#########
# Const #
#########
declare -r DEFAULT_PRIMER='ATCTCTCTCAATTTTTTTTTTTTTTTTTTTTTTTAAGAGAGAGAT'
declare -ri DEFAULT_NUM_THREADS=8

#########
# Usage #
#########
function usage {
cat << EOF
${PACKAGE_NAME}::${MODULE_NAME}
================

Split PacBio subreads BAM file with a custom smrtbell sequence

OPTIONS:
        -h      Show usage
        -v      Print version

${REQUIRED}[ required ]
        -i      input files with pacBio canonical smrtbell adaptor clipped, can be either of the following two types.
                1. [ RSII ] directory to RSII primary analysis output, such as A01_1
                2. [ Sequel ] subreads.bam
                You can specify multiple inputs, each with its own -i tag.
                Avoid using inputs with the same basename.

${OPTIONAL}[ optional ]
        -p      the custom smrtbell primer sequence, Default: ${DEFAULT_PRIMER}
        -o      Output folder, will create if not exist. Default: ${PWD}
        -t      Number of threads to use. Default: ${DEFAULT_NUM_THREADS}

${ADVANCED}[ advanced ]
        -D      use debug mode (bash -x). Default: off

EOF
echo -e "${FONT_COLOR_RESET}"
}

##########
# Config #
##########
declare -a REQUIRED_PROGRAMS=('bax2bam' 'split_primer_from_pbbam')
declare -a Inputs=()

#############################
# ARGS reading and checking #
#############################
while getopts "hvi:p:o:t:D" OPTION; do
    case $OPTION in
        h)  usage && exit 0 ;;
        v)  echo ${PACKAGE_NAME}::${MODULE_NAME} v${MODULE_VERSION} && exit 0;;
        i)  declare Input=$(readlink -f ${OPTARG});
            Inputs+=("${Input}")
            ;;
        p)  declare PrimerSequence=${OPTARG};;
        o)  declare OutputDir=${OPTARG};;
        t)  declare -i Threads=${OPTARG};;
        D)  set -x;;
        *)  usage && exit 1;;
    esac
done

[[ -z $Input ]] && usage && echo2 "Need do provide at least one input fl_nc fasta file with -i option" error
[[ -z $PrimerSequence ]] && declare PrimerSequence=$DEFAULT_PRIMER && echo2 "You have not provide the pipeline with a custom smrtbell primer with -p, thus to use the default primer ${DEFAULT_PRIMER}" warning
[[ -z $Threads ]] && declare -i Threads=${DEFAULT_NUM_THREADS}
[[ -z $OutputDir ]] && OutputDir=${RANDOM}.out

for file in "${Inputs[@]}"; do dirOrFileCheck $file; done
for program in "${REQUIRED_PROGRAMS[@]}"; do binCheck $program; done

mkdir -p $OutputDir || echo2 "Don't have permission to create folder $OutputDir" error
cd $OutputDir || echo2 "Don't have permission to access folder $OutputDir" error
( touch a && rm -f a ) || echo2 "Don't have permission to write in folder $OutputDir" error

declare RUNUID=`echo $@ | md5sum | cut -d' ' -f1` # for re-run
echo2 "Begin the pipeline ${PACKAGE_NAME} ${MODULE_NAME}"
declare -i Step=1

echo2 "Validate the input files"
declare -a bamInputs

for input in "${Inputs[@]}"; do
    declare Prefix=$(basename ${input})
    echo2 "Validating input ${input}"
    if isDir ${input}; then
        echo2 "Input is an directory, check for RS II data"
        if ! validateRSII $input; then echo2 "invalid RS II output $input" error; fi
        echo2 "Converting legacy RS II data to bam format"
        declare Bax2BamArgs=$(find ${input}/Analysis_Results -name "*bax.h5"  | xargs echo)
        if ! isFile ${Prefix}.subreads.bam; then
            bax2bam $Bax2BamArgs -o ${Prefix} || echo2 "failed to convert bax back to bam format" error
        else
            echo2 "Skipping ${input} because the output ${Prefix}.subreads.bam has existed" warning
        fi
        bamInputs+=( "${Prefix}.subreads.bam" )
    else
        bamInputs+=( "${input}" )
    fi
done
let Step+=1

echo2 "Begin to split bam files"
if [[ ! -f ${RUNUID}.${Step}.Done || ${RUNUID}.${Step}.Done -ot ${RUNUID}.$((Step-1)).Done ]]; then
    for bamFile in "${bamInputs[@]}"; do
        declare Prefix=$(basename ${bamFile%.subreads.bam})
        echo2 "Split ${bamFile}"
        if ! isFile ${Prefix}.refarm.bam; then
            split_primer_from_pbbam -p $PrimerSequence -o ${Prefix}.refarm.bam -t ${Threads} ${bamFile}
        else
            echo2 "Skipping ${bamFile} because the output ${Prefix}.refarm.bam has existed" warning
        fi
    done
fi
let Step+=1
