#!/bin/bash
#set -o xtrace
date_now=`date +%Y-%m-%d.%H.%M`
new_database_name=gp_`date +%Y%m%d`
db_name=gp_dev
db_con=METASCAPE
email_addr=xxx@xxx.org
config=config.xml
log=logs/build_dev_log_$date_now.txt

if true; then
    mkdir -p upload_backup
    mkdir -p download_backup
    mv download download_backup/download_$date_now
    mv upload upload_backup/upload_$date_now
    mkdir -p download
    mkdir -p upload
    mkdir -p logs
    mkdir -p logs/error_logs
fi

#backup the last statistics files
if true;then
    folder_keep=~/keep_$date_now
    mkdir $folder_keep
    cp /var/www/html/gp/Content/MenuPages/statistics.html   $folder_keep
    cp -r /var/www/html/gp/Content/MenuPages/out   $folder_keep
    cp /var/www/html/gp/Content/MenuPages/gp_description_stats.js   $folder_keep
fi


#main branch, sync all the data source
if true; then
    cmd="python -u gp.py --rebuild N --skipvalidation -c $config -n $db_con -d $db_name -rf download_backup/download_$date_now -el logs/error_logs/gp_dev_$data_now.txt > $log"
    printf "$cmd\n"
    time eval "$cmd"
fi

if true; then
    cmd="python -u gp.py --history  --connection $db_con --database $db_name >> $log"
    printf "$cmd\n"
    eval "$cmd"

    cmd="python -u gp.py --makestatisticshtml --connection $db_con --database $db_name >> $log"
    printf "$cmd\n"
    eval "$cmd"
fi

if true;then
    fn="~/sync/$new_database_name.sql"
    cmd="./db_backup.py -c $db_con -d $db_name -o $fn -a backup"
    printf "$cmd\n"
    eval "$cmd"
    cmd="./db_backup.py -i $fn -r $new_database_name -c $db_con -a restore"
    printf "$cmd\n"
    eval "$cmd"
fi

