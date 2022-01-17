#!/bin/bash
git log --all --format='%aN <%cE>' | sort -u > AUTHORS
