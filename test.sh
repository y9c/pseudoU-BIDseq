#!/usr/bin/env bash
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2022-04-04 18:49

logo="
\033[1;34m  ██╗ ██╗  ██╗       ██████╗ ██╗██████╗ ███████╗███████╗ ██████╗  \033[0m
\033[1;34m  ██║ ██║  ██║       ██╔══██╗██║██╔══██╗██╔════╝██╔════╝██╔═══██╗ \033[0m
\033[1;34m  ╚█████║███╔╝ █████╗██████╔╝██║██║  ██║███████╗█████╗  ██║   ██║ \033[0m
\033[1;34m   ╚══██║═══╝  ╚════╝██╔══██╗██║██║  ██║╚════██║██╔══╝  ██║▄▄ ██║ \033[0m
\033[1;34m      ██║            ██████╔╝██║██████╔╝███████║███████╗╚██████╔╝ \033[0m
\033[1;34m      ╚═╝            ╚═════╝ ╚═╝╚═════╝ ╚══════╝╚══════╝ ╚══▀▀═╝  \033[0m
\033[1;34m                                                                  \033[0m
\033[1;34m      ██████╗ ██╗██████╗ ███████╗██╗     ██╗███╗   ██╗███████╗    \033[0m
\033[1;34m      ██╔══██╗██║██╔══██╗██╔════╝██║     ██║████╗  ██║██╔════╝    \033[0m
\033[1;34m      ██████╔╝██║██████╔╝█████╗  ██║     ██║██╔██╗ ██║█████╗      \033[0m
\033[1;34m      ██╔═══╝ ██║██╔═══╝ ██╔══╝  ██║     ██║██║╚██╗██║██╔══╝      \033[0m
\033[1;34m      ██║     ██║██║     ███████╗███████╗██║██║ ╚████║███████╗    \033[0m
\033[1;34m      ╚═╝     ╚═╝╚═╝     ╚══════╝╚══════╝╚═╝╚═╝  ╚═══╝╚══════╝    \033[0m
"

printf "$logo\n"

# Set default values
snake="/opt/pipeline/Snakefile"
conf="data.yaml"
cores=$(python -c 'import yaml;print(yaml.safe_load(open("data.yaml","r")).get("cores",""))')
if [ -z "$cores" ]; then cores=48; fi

POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
  case $1 in
  -j | --jobs | --cores)
    EXTENSION="$2"
    shift # past argument
    shift # past value
    ;;
  -c | --conf)
    conf="$2"
    shift # past argument
    shift # past value
    ;;
  -s | --snake)
    snake="$2"
    shift # past argument
    shift # past value
    ;;
  *)
    POSITIONAL_ARGS+=("$1") # save positional arg
    shift                   # past argument
    ;;
  esac
done

echo $@
