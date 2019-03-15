#!/usr/bin/awk -f

BEGIN{
    flag=1
    header=1
}

/^$|(HepMC::(Version|IO_GenEvent))/{
    flag=0
}

{
    if (flag==1 || header==1) print
}

/HepMC::IO_GenEvent-START_EVENT_LISTING/{
    header=0
    flag=1
}

/HepMC::IO_GenEvent-END_EVENT_LISTING/{
    footer=$0
}

END{
    print footer
}
