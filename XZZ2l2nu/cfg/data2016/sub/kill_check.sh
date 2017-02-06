#!/bin/sh


kill `ps aux | grep check.sh | awk {'print $2'}`

