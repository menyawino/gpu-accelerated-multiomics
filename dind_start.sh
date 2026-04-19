#!/bin/bash
dockerd --iptables=false --ip-forward=false --bridge=none &> /tmp/dockerd.log &
