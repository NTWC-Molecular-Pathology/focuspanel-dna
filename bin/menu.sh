#!/bin/bash

trap ' ' 2 15 20
python menu.py
trap - 2 15 20