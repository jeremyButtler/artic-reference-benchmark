 { # Main
   if($4 ~ /[Aa][Ll][Tt]/) next; # Alternate primer
   if($4 ~ /[Ll][Ee][Ff][Tt]/)
   { # If have a left primer
     ++numPrimI;
     lPrimAry[numPrimI] = $4;
     sub(/[^_]*./, "", $4); # remove prefix
     sub(/_.*/, "", $4);    # remove ending
     lprimNumAry[numPrimI] = $4; # primer number
   } # If have a left primer

   else
   { # else if this is a right primer
      for(iPrim = 1; iPrim <= numPrimI; ++iPrim)
      { # Loop: Look for matching primer
         primNumI = $4;
         sub(/[^_]*./, "", primNumI); # remove prefix
         sub(/_.*/, "", primNumI);    # remove ending

         if(lprimNumAry[iPrim] == primNumI)
         { # If I found the match
            printf "%s\t%s\n", lPrimAry[iPrim], $4;
            next; # Restart on the next line
         } # If I found the match
      } # Loop: Look for matching primer
   } # else if this is a right primer
} # MAIN

