#!/bin/bash

LOG_FILE="memory_usage_deming.log"
check_memory() {
    total_memory=$(free | awk '/^Mem:/ {print $2}')
    used_memory=$(free | awk '/^Mem:/ {print $3}')
    
    usage_percent=$(( (used_memory * 100) / total_memory ))
    
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    echo "$timestamp - Memory Usage: $usage_percent%" >> "$LOG_FILE"
    
    if [ "$usage_percent" -gt 90 ]; then
        echo "$timestamp - Memory usage exceeds 95%. Killing all screen sessions..." >> "$LOG_FILE"
        pkill screen
    fi
}

while true; do
    check_memory
    sleep 60
done