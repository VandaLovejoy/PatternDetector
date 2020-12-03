import java.util.*; import java.io.*; import java.math.* ; import java.lang.*;;

public class ScanItFast implements Runnable {
    static boolean 	VERBOSE = false ,
            PRINTALL = false ;
    private boolean [] hasChars , keepMe, isAmbiguous, isNotUnique ;
    private String[] mafTab ;
    private String[] key;
    private ArrayList <String[]> motifs;
    ArrayList<char[]> alnTab;
   // private char[][] AlnTab;
    private static String Path, SSZBINARY, ALIFOLDBINARY;
    private int [] coordTab;
    private int  goodSeqs, iterate, GAPS, retainedColumns;
    private BufferedWriter WriteALN ;
    private BufferedReader ReadFile;
    private String [] OutAln, OutAlnRC, FilteredTab,	NameTab,  TempTab = new String [1] ;
    private String Line = "" ;
    private File Aln, AlnRC ;
    private double [] stats, chars, totalChars ;
    private static double	mpi = 0,
            mpi2 = 0,
            var = 0,
            shanon = 0,
            uniqueComps = 0,
            uniqueSeqs,
            SSZ_THRESHOLD = -2.7 , 		// alignments scoring below this will be kept (Z-score)
            SSZR_THRESHOLD = -2.2 ,		// alignments scoring below this will be kept (Z-score)
            RNAZ_THRESHOLD = 0.32 ,   	// alignments scoring above this will be kept (SVM probability)
            outCols ;
    private double [][]	pids, gaps ;

    ScanItFast(String[] mafTab, ArrayList motifs, ArrayList<char[]> alnTab,
               String[] key, String Path, int GAPS,
               String SSZBINARY,  boolean VERBOSE, boolean PRINTALL ) {
        this.mafTab = mafTab ;
        this.Path = Path ;
        this.SSZBINARY = SSZBINARY ;
        this.VERBOSE = VERBOSE ;
        this.PRINTALL = PRINTALL ;
        this.alnTab = alnTab ;
        this.GAPS = GAPS ;
        this.motifs = motifs;
        this.alnTab = alnTab;
        this.key = key;
    }

    public void run() {
        if (VERBOSE)
            System.out.println("- - -> Starting Scan") ;
        if (VERBOSE && alnTab.size() != motifs.size() ) {
            System.out.println(" #### Maf and AlnTab aren't same length" ) ;
            System.out.println( motifs.size() +" "+ alnTab.size() ) ;
        }
        isAmbiguous = new boolean [ alnTab.size() ];
        int startPos = Integer.parseInt(mafTab[2]);
        // remove identical rows or those with too many gaps & N's
        Set<String>	Uniques = new LinkedHashSet<String>(),
                UniquesWithGaps = new LinkedHashSet<String>(),
                UniqueNames = new LinkedHashSet<String>() ;
        for (int seq = 0 ; seq != alnTab.size() ; seq++ ) {
            String DegappedWindow = new String( alnTab.get(seq)).replaceAll("[^ATCGUatcgu]", "" ).toUpperCase() ;
            // only retains non-identical unaligned sequences with at least one character
            if ( DegappedWindow.length() > 0 && Uniques.add( DegappedWindow ) ) {
                UniquesWithGaps.add(  new String( alnTab.get(seq)).toUpperCase() ) ;
                UniqueNames.add(  motifs.get(seq)[ 0 ] ) ;
            }
        }
            FilteredTab = new String [ UniquesWithGaps.size()]; // creating an Array from Hash
            FilteredTab = UniquesWithGaps.toArray(FilteredTab);
            NameTab = new String [ UniquesWithGaps.size() ] ;
            NameTab = UniqueNames.toArray( NameTab ) ;
// first check for > 2 seqs
        goodSeqs = UniquesWithGaps.size()  ;
        if ( goodSeqs <= 3 ) {
            if (VERBOSE)
                System.out.println("-> Not Enough seqs in this window!") ;
            return;
        }

        // remove gappy sequences
        if (VERBOSE)
            System.out.println("- -> Gappy sequences") ;
        keepMe = new boolean [ FilteredTab.length ] ;
        for (int seq = 0 ; seq != FilteredTab.length ; seq++ ) {
           /*System.out.println("this is filtered tab with gaps:" + FilteredTab[seq].length() + "\n" +
           "this is filtered tab without gaps"+  FilteredTab[ seq ].replaceAll("[^ATCGUatcgu]", "" ).length()+
           "\n" + " this is the limit that has to be reached " + (int)(FilteredTab[seq].length()*((double)50/100)));*/
            if ( FilteredTab[ seq ].replaceAll("[^ATCGUatcgu]", "" ).length() >= (int)(FilteredTab[seq].length()*((double)50/100))) {
                keepMe[seq] = true ;
            }
            else {
                keepMe[seq] = false ;
                goodSeqs-- ;
                //if (VERBOSE)
                //	System.out.println("  --> removed a GAPpy sequence form the alignment" ) ;
            }

        }
       // System.out.println("this is the end of block");
        if ( goodSeqs <= 3 ) {
            if (VERBOSE)
                System.out.println("-> Not Enough seqs in this window!") ;
            return;
        }
        // exit when human is shit
        if ( !keepMe[ 0 ] )
            return;

        // check for gap-only columns
        if (VERBOSE)
            System.out.println("- -> Gap only columns") ;
        retainedColumns = FilteredTab[0].length();
        hasChars = new boolean [ FilteredTab[0].length() ];
        gapScan: for (int col = 0 ; col < FilteredTab[0].length() - 1 ; col++ ) {
            for (int seq = 0 ; seq != FilteredTab.length ; seq++ ) {
                if ( keepMe[ seq ] ) {
                    if ( FilteredTab[ seq ].charAt( col ) == 'A' || FilteredTab[ seq ].charAt( col ) == 'C'
                            || FilteredTab[ seq ].charAt( col ) == 'T' || FilteredTab[ seq ].charAt( col ) == 'G' ) {
                        hasChars[ col ] = true ;
                        continue gapScan ;
                    }
                }
            }
            if ( !hasChars[ col ] )
                retainedColumns-- ;
            //if ( !hasChars[ col ] && VERBOSE)
            //	System.out.println( "-> empty col!" );
        }

        // prepare clustalw file
        if (VERBOSE)
            System.out.println("- -> preparing Clustal format") ;
        OutAln = new String[ goodSeqs ];
        OutAlnRC = new String[ goodSeqs ]  ;
        iterate = 0 ;
        for (int seq = 0 ; seq != FilteredTab.length  ; seq++ ) { //removed x < goodseqs
            if ( keepMe[ seq ] ) {
                OutAln[iterate] = NameTab[seq].substring(0, Math.min(NameTab[seq].length(), 20));
                for (int i = 0; i != 25 - Math.min(NameTab[seq].length(), 20); i++)
                    OutAln[iterate] = OutAln[iterate] + " ";
                for (int i = 0; i != FilteredTab[0].length(); i++)
                    if (hasChars[i])
                        OutAln[iterate] = OutAln[iterate] + FilteredTab[seq].charAt(i);
                OutAln[iterate] = OutAln[iterate] + "\n";
                iterate++;

            }
        }
//*********************************************************************
        //					calculate stats						*
        //*********************************************************************
        if (VERBOSE)
            System.out.println("- - -> calculating statistics") ;
        uniqueSeqs = goodSeqs;
        outCols = OutAln[0].length()-25 ; //change last variable if CLUSTAL properties changes
        stats = new double [6];
        chars = new double [5];
        totalChars = new double [5];
        pids = new double [ goodSeqs ][ goodSeqs ];
        gaps = new double [ goodSeqs ][ goodSeqs ]; // gaps (and potentially mismatches)
        isNotUnique = new boolean [ goodSeqs ] ;
        // calculate id matrix and mpi
        for ( int k = 25 ; k != OutAln[0].length()-1 ; k ++ ) { // -1 avoids line break char, 25 is clustal seq ID length
            lines :for ( int i = 0 ; i != goodSeqs ; i++ ) {
                // initiate gaps[] and pids[] to 0 ???????????????????????
                if ( isNotUnique[ i ] )
                    continue lines ;
                for ( int j = i+1 ; j != goodSeqs ; j++ ) {
                    try {

                        if ( isNotUnique[ j ] )
                            continue;
                        else if ( OutAln[ i ].charAt( k ) == OutAln[ j ].charAt( k ) ) {
                            // this DP matrix makes shit easy!
                            if ( OutAln[ i ].charAt( k ) == 'A' || OutAln[ i ].charAt( k ) == 'T'
                                    || OutAln[ i ].charAt( k ) == 'C' || OutAln[ i ].charAt( k ) == 'G'
                                    || OutAln[ i ].charAt( k ) == 'U') { // U is just in case alignments are RNA
                                pids[ i ][ j ]++ ;
                                pids[ j ][ i ]++ ;
                            }
                        }
                        // this ignores "-:-"
                        else if ( OutAln[ i ].charAt( k ) != OutAln[ j ].charAt( k ) ) {
                            pids[ j ][ i ]++ ; // mismatch
                            if ( OutAln[ i ].charAt( k ) == '-' || OutAln[ i ].charAt( k ) == '.')
                                gaps[i][j]++; // gap
                            if ( OutAln[ j ].charAt( k ) == '-' || OutAln[ j ].charAt( k ) == '.')
                                gaps[j][i]++; // gap
                        }
                    }catch (Exception E) {
                        E.printStackTrace();
                        System.err.println( "Caught Exception!\n");
                        System.err.println( i +" " +j +" " + k + " " + OutAln.length + " " + goodSeqs
                                + " FilteredTab[i]=" + FilteredTab[i].length() + " FilteredTab[j]=" + FilteredTab[j].length());
                        System.err.println( "OutAln[i]=" + OutAln[i].length() + " OutAln[j]=" + OutAln[j].length()
                                + "\n" + OutAln[i] + "\n" + OutAln[j]) ;
                    }
                    // keep unique seqs ignoring gaps
                    if ( k == OutAln[0].length()-2 ) {
                        if (	pids[ j ][ i ]  - gaps[ i ][ j ] ==  pids[ i ][ j ]  ||
                                pids[ j ][ i ]  - gaps[ j ][ j ] ==  pids[ i ][ j ]  ) {
                            //both sequences are identical without gaps
                            //keep the longer one
                            if ( gaps[ i ][ j ] > gaps[ j ][ i ] )
                                isNotUnique[ i ] = true ;
                            else
                                isNotUnique[ j ] = true ; // this should also consider identical seqs
                        }
                        else {
                            uniqueComps++ ;
                            // old mean pairwise identity ( considers gaps )
                            mpi = mpi + 100 * pids[ i ][ j ] / pids[ j ][ i ] ;
                            // classical average identity
                            mpi2 = mpi2 + 100 * pids[ i ][ j ] / Math.min( OutAln[ i ].replaceAll("[^ATCGU]","").length(),
                                    OutAln[ j ].replaceAll("[^ATCGU]","").length() ) ;
                        }
                    }
                }
            }
        }
        // calculate gaps, GC, shanon entropy, and Reverse Complement
        for ( int k = 25 ; k != OutAln[0].length() ; k ++ ) {
            chars = new double [5] ;
            for ( int i = 0 ; i != goodSeqs ; i++ ) {
                if (isNotUnique[i]) {
                    if (k == OutAln[0].length() - 2)
                        uniqueSeqs--;
                    continue;
                }

                switch (OutAln[i].charAt(k)) {
                    case 'A':
                        chars[0]++;
                        totalChars[0]++;
                        OutAlnRC[i] = (k == 25) ? "T" : "T" + OutAlnRC[i];
                        break;
                    case 'U':
                        chars[1]++;
                        totalChars[1]++;
                        OutAlnRC[i] = (k == 25) ? "A" : "A" + OutAlnRC[i];
                        break;
                    case 'T':
                        chars[1]++;
                        totalChars[1]++;
                        OutAlnRC[i] = (k == 25) ? "A" : "A" + OutAlnRC[i];
                        break;
                    case 'C':
                        chars[2]++;
                        totalChars[2]++;
                        OutAlnRC[i] = (k == 25) ? "G" : "G" + OutAlnRC[i];
                        break;
                    case 'G':
                        chars[3]++;
                        totalChars[3]++;
                        OutAlnRC[i] = (k == 25) ? "C" : "C" + OutAlnRC[i];
                        break;
                    case '\n':
                        OutAlnRC[i] = OutAlnRC[i] + '\n';
                        break;
                    case 'N':
                        chars[4]++;
                        totalChars[4]++;
                        OutAlnRC[i] = (k == 25) ? "N" : "N" + OutAlnRC[i];
                        break;
                    default:
                        chars[4]++;
                        totalChars[4]++;
                        OutAlnRC[i] = (k == 25) ? "-" : "-" + OutAlnRC[i];
                        break;
                }
            }
            for (int z = 0 ; z != 5 ; z++ )
                shanon = ( chars[z] == 0 )? shanon + 0 : shanon +  chars[z]/uniqueSeqs * ( Math.log( chars[z]/uniqueSeqs ) / Math.log( 2 ));
        }
        //System.out.println( uniqueSeqs +"\t"+goodSeqs+"\t"+outCols+"\t"+totalChars[4]+"\t"+( outCols * goodSeqs));
        stats[0] = mpi / uniqueComps;																		// Mean Pairwise ID
        stats[5] = mpi2 / uniqueComps;																    // classical MPI
        for (int seq1 = 0 ; seq1 != goodSeqs ; seq1++ )
            for (int seq2 = seq1 +1 ; seq2 != goodSeqs ; seq2++ )
                if ( !isNotUnique[ seq1 ] && !isNotUnique[ seq2 ] )
                    var = var + (double) Math.pow( (100*pids[ seq1 ][ seq2 ]/ pids[ seq2 ][ seq1 ]) - stats[0] , 2) ;
        stats[1] = var / uniqueComps ;																// Variance
        stats[2] = -1 * shanon / ((double)outCols) ;													    // Normalized Shanon entropy
        stats[3] = 100*(totalChars[2]+totalChars[3])/(totalChars[0]+totalChars[1]+totalChars[2]+totalChars[3]) ;	   // GC content
        stats[4] = 100 * totalChars[4] / (outCols * uniqueSeqs) ;										  // GAP content
       // System.out.println( stats[0]+"\t"+(Math.sqrt(stats[1]))+"\t"+stats[2]+"\t"+stats[3]+"\t"+stats[4]) ; 	 // print stats
        System.out.println( stats[0]);
        // save BED coords from MAF file
        if (VERBOSE)
            System.out.println("- -> Calculating BED coords ") ;
        String	BedFile = mafTab[1].substring( mafTab[1].lastIndexOf(".")+1) +"\t";
        if ( mafTab[4].equals("+") ){
            BedFile = BedFile + (startPos) + Integer.parseInt(key[0]) +"\t"+(startPos + Integer.parseInt(key[1]))+"\t" ;
        }
        else { // this should only occur in user specified cases
            BedFile = BedFile + (Integer.parseInt(mafTab[5]) - (startPos + Integer.parseInt(key[0])) - Integer.parseInt(key[1])
                    +"\t"+ (Integer.parseInt(mafTab[5]) - (startPos+ Integer.parseInt(key[0]) +1 )))+"\t";
        }
        BedFile = BedFile	+ (int)uniqueSeqs+":"+((double)(int)(10*stats[0])/10)+":"      // MPI
                + ((double)(int)(10*stats[5])/10)+":"					  // CLASSIC MPI
                + ((double)(int)(10*stats[4])/10) +":"					 // GAPS
                + ((double)(int)(10*Math.sqrt(stats[1]))/10) +":"			// STDEV
                + ((double)(int)(100*stats[2])/100)  +":"			    // SHANON
                + ((double)(int)(10*stats[3])/10) ;				   //      GC
        if (VERBOSE)
            System.out.println( "Pre SISSIz bed file: \n"+" "+BedFile ) ;

        int random = (int)((double)10000*Math.random()) ;
        File Aln = new File( Path+"/"+BedFile.replaceAll("\t","_")+".aln."+ random ),	//
                AlnRC = new File( Path+"/"+BedFile.replaceAll("\t","_")+"rc.aln." + random );  //
        // v v v v v v v v    INCLUSION STATS     v v v v v v v v v v v v v
        if ( outCols > (FilteredTab[0].length()) / 2 && uniqueSeqs > 2 &&  stats[4] <= 75 && stats[0] > 0 ) {
            // Write Sequences to ALN Format
            /*try {
                BufferedWriter WriteClustal = new BufferedWriter(new FileWriter( Aln )),
                        WriteClustalRC = new BufferedWriter(new FileWriter( AlnRC ));
                WriteClustal.write("CLUSTAL format sucks\n\n") ;
                WriteClustalRC.write("CLUSTAL format sucks\n\n") ;
                for ( int y = 0 ; y != goodSeqs ; y++ ) {
                    if ( !isNotUnique[ y ] ) {
                        WriteClustal.write( OutAln[ y ] ) ;
                        OutAlnRC[ y ] = OutAln[ y ].substring(0,25)+ OutAlnRC[ y ]  ;
                        WriteClustalRC.write( OutAlnRC[ y ] ) ;
                    }
                }
                WriteClustal.close() ;
                WriteClustalRC.close() ;
            } catch (IOException Fuck) {
                if (VERBOSE)
                    System.err.println("Arrgh... Couldn't write clustal file!");
                Fuck.printStackTrace();
                Aln.delete() ;
                AlnRC.delete() ;
                return;
            }*/
        }
        else {
            if (VERBOSE) {
                System.out.println("---> rejected alignment" ) ;
                System.out.println("     outcols = "+outCols +"\tuniqueseqs = "+uniqueSeqs+
                        "\tGAPS = "+stats[4]+"\n    PID = "+stats[0]);
                if ( stats[0] < 5 )
                    System.out.println("-----> SUPER LOW PID");
            }
			Aln.delete() ;
			AlnRC.delete() ;
            return ;
        }
        String	FinalBedFile = "",
                FinalBedFileRC = "",
                Antisense = (mafTab[4].equals("+"))? "-" : "+" ;
//***************** 	RNALalifold scan & parse		******************
        /*if ( stats[0] <= 85 ) { // 85% sequence identity cut-off for RNAz
            String[] SissizOutTab = new String[12];
            try {
                SissizOutTab = ScanSSZ(Path, BedFile, stats, 1, random);
                if (SissizOutTab == null) { // timeout
                    Aln.delete();
                }
            } catch (IOException Fuck) {
                Fuck.printStackTrace();
                System.err.println("ScanSSZ failed with ");
                for (int y = 0; y != goodSeqs; y++) {
                    if (!isNotUnique[y]) {
                        System.err.println(OutAln[y]);
                    }
                }
                Aln.delete();
            }
        }*/
/*System.out.println("this is the mafTab: " + mafTab[2] +"\n" +

        "this is the interval of motif" + Arrays.toString(key) + "\n"+ "this is the number of aln:" + motifs.size() );*/

    }
}
