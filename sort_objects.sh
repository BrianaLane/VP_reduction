#!/bin/bash

SCI="object"
morn="morn eve"
#morn="morn"
#OBJKEYS="M82_sky M82_F1 M82_F2 Faige34 NGC4666_sky NGC4666 Feige67" #VP_data
OBJKEYS="M82_Fa M82_F7 M82_F8 M82_F9 M82_F10 M82_F11 M82_F12 M82_F13 M82_F14 M82_sky2 GRW70d5824 Feige34" #VP_data2
#OBJKEYS="M82_sky M82_F1 M82_F2 Feige34 NGC3079 NGC3079_sky NGC4631_sky NGC4631 Feige67" #VP_data3
#OBJKEYS="Feige34 NGC3079_sky NGC3079 Feige67" #VP_data4

FT=$CUREVP

getkeyword()
{
    # extract value of keyword $2 from file $1
    #$FT/headfits $1 | grep $2 | cut -b 12-19 | sed 's/ //g'
    $CUREVP/headfits -r -k $2 $1
}


# Check that the necessary environment variables are set
if [ x$CUREVP == "x" ]; then
    echo "Please set CUREVP environment variable. e.g."
    echo "bash : export CUREVP=~/cure/bin
    echo "csh  : setenv CUREVP "~/cure/bin"
        exit 1
fi

# Check that subtractfits is in place (we will assume that the rest
#  of the fitstools is as well if subtractfits is.
if [ ! -x $CUREVP/headfits ]; then
    echo "Unable to find $CUREVP/headfits, please check that CUREVP is"
    echo "  set correctly."
    exit 1
fi


for month in $@; do

    for day in `cd redux/$month; echo 20??????`; do

        mkdir -p data/flat/$day
        mkdir -p data/trace/$day

        for d in $morn $morn; do

            for objframe in `ls redux/$month/$day/object/$d/Spes*.fits`; do
                
                objname=`getkeyword  $objframe "OBJECT" | tr -s ' ' | cut -d ' ' -f 1`
                file=`basename $objframe .fits | sed -e 's/S//'`

                cont=0
                for ob in $OBJKEYS; do
                    echo $objname | grep $ob &> /dev/null
                    if [ $? -eq 0 ]; then
                        cont=1
                        break
                    fi
                done

                if [ $cont -eq 0 ]; then 
                    echo "Ignoring $objname"
                    continue
                fi

#                echo "Still here $objname"
                 
                mkdir -p data/$objname
                
                cp redux/$month/$day/object/$d/S${file}.fits data/$objname/
                cp redux/$month/$day/object/$d/e.S${file}.fits data/$objname/
                
                #cp redux/$month/$day/flat/$d/n$file.fits data/flat/$day/
                cp redux/$month/$day/flat/$d/$file.fmod data/flat/$day/
                ln -fs ../flat/$day/$file.fmod data/$objname/
                cp redux/$month/$day/flat/$d/$file.dist data/flat/$day/
                ln -fs ../flat/$day/$file.dist data/$objname/
                
                cp redux/$month/$day/${d}_mastertrace.pmod data/trace/$day/${d}_mastertrace.pmod
                
                ln -fs ../trace/$day/${d}_mastertrace.pmod data/$objname/$file.pmod

            done

        done

    done
done
