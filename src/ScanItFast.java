import java.util.*; import java.io.*; import java.math.* ; import java.lang.*;;

public class ScanItFast implements Runnable {
    static boolean VERBOSE = false,
            PRINTALL = false;
    private boolean[] hasChars, keepMe, isAmbiguous, isNotUnique;
    private String[] key;
    private ArrayList<String[]> motifs;
    ArrayList<char[]> alnTab;
    private static String Path, dirPath, SSZBINARY, RSCAPEBINARY;
    private int[] coordTab;
    private int goodSeqs, iterate, GAPS, retainedColumns;
    private BufferedWriter WriteALN;
    private BufferedReader ReadFile;
    private String[] OutAln, OutAlnRC, FilteredTab, NameTab, TempTab = new String[1];
    private String Line = "";
    private File Aln, AlnRC;
    private double[] stats, chars, totalChars;
    private static double mpi = 0,
            mpi2 = 0,
            var = 0,
            shanon = 0,
            uniqueComps = 0,
            uniqueSeqs,
            SSZR_THRESHOLD = -3.0,        // alignments scoring below this will be kept (Z-score)
            outCols;
    private double[][] pids, gaps;

    ScanItFast(ArrayList motifs, ArrayList<char[]> alnTab,
               String[] key, String Path, String dirPath, int GAPS,
               String SSZBINARY, String RSCAPEBINARY, boolean VERBOSE, boolean PRINTALL) {
        this.Path = Path;
        this.dirPath = dirPath;
        this.SSZBINARY = SSZBINARY;
        this.RSCAPEBINARY = RSCAPEBINARY;
        this.VERBOSE = VERBOSE;
        this.PRINTALL = PRINTALL;
        this.alnTab = alnTab;
        this.GAPS = GAPS;
        this.motifs = motifs;
        this.alnTab = alnTab;
        this.key = key;
    }

    public void run() {
        if (VERBOSE)
            System.out.println("- - -> Starting Scan");
        if (VERBOSE && alnTab.size() != motifs.size()) {
            System.out.println(" #### Maf and AlnTab aren't same length");
            System.out.println(motifs.size() + " " + alnTab.size());
        }
        isAmbiguous = new boolean[alnTab.size()];
        int startPos = Integer.parseInt(key[1]);
        // remove identical rows or those with too many gaps & N's

        ArrayList<String> Uniques = new ArrayList();
        ArrayList<String> UniquesWithGaps = new ArrayList();
        Set<String> UniqueNames = new LinkedHashSet<String>();
        for (int seq = 0; seq != alnTab.size(); seq++) {
            String DegappedWindow = new String(alnTab.get(seq)).replaceAll("[^ATCGUatcgu]", "").toUpperCase();
            // only retains non-identical unaligned sequences with at least one character
            if (DegappedWindow.length() > 0 && UniqueNames.add(motifs.get(seq)[0])) {
                Uniques.add(DegappedWindow);
                UniquesWithGaps.add(new String(alnTab.get(seq)).toUpperCase());
            }
        }

        FilteredTab = new String[UniquesWithGaps.size()]; // creating an Array from Hash
        FilteredTab = UniquesWithGaps.toArray(FilteredTab);
        NameTab = new String[UniquesWithGaps.size()];
        NameTab = UniqueNames.toArray(NameTab);
// first check for > 2 seqs
        goodSeqs = UniquesWithGaps.size();
        if (goodSeqs <= 3) {
            if (VERBOSE)
                System.out.println("-> Not Enough seqs in this window!");
            return;
        }

        // remove gappy sequences
        if (VERBOSE)
            System.out.println("- -> Gappy sequences");
        keepMe = new boolean[FilteredTab.length];
        for (int seq = 0; seq != FilteredTab.length; seq++) {
           /*System.out.println("this is filtered tab with gaps:" + FilteredTab[seq].length() + "\n" +
           "this is filtered tab without gaps"+  FilteredTab[ seq ].replaceAll("[^ATCGUatcgu]", "" ).length()+
           "\n" + " this is the limit that has to be reached " + (int)(FilteredTab[seq].length()*((double)50/100)));*/
            if (FilteredTab[seq].replaceAll("[^ATCGUatcgu]", "").length() >= (int) (FilteredTab[seq].length() * ((double) 50 / 100))) {
                keepMe[seq] = true;
            } else {
                keepMe[seq] = false;
                goodSeqs--;
            }

        }
        // System.out.println("this is the end of block");
        if (goodSeqs <= 3) {
            if (VERBOSE)
                System.out.println("-> Not Enough seqs in this window!");
            return;
        }
        // exit when human is shit
        if (!keepMe[0])
            return;

        // check for gap-only columns
        if (VERBOSE)
            System.out.println("- -> Gap only columns");
        retainedColumns = FilteredTab[0].length();
        hasChars = new boolean[FilteredTab[0].length()];
        gapScan:
        for (int col = 0; col < FilteredTab[0].length(); col++) {
            for (int seq = 0; seq != FilteredTab.length; seq++) {
                if (keepMe[seq]) {
                    if (FilteredTab[seq].charAt(col) == 'A' || FilteredTab[seq].charAt(col) == 'C'
                            || FilteredTab[seq].charAt(col) == 'T' || FilteredTab[seq].charAt(col) == 'G') {
                        hasChars[col] = true;
                        continue gapScan;
                    }
                }
            }
            if (!hasChars[col])
                retainedColumns--;
        }

        // prepare clustalw file
        if (VERBOSE)
            System.out.println("- -> preparing Clustal format");
        OutAln = new String[goodSeqs];
        OutAlnRC = new String[goodSeqs];
        iterate = 0;
        for (int seq = 0; seq != FilteredTab.length; seq++) { //removed x < goodseqs
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
            System.out.println("- - -> calculating statistics");
        uniqueSeqs = goodSeqs;
        outCols = OutAln[0].length() - 25; //change last variable if CLUSTAL properties changes
        stats = new double[6];
        chars = new double[5];
        totalChars = new double[5];
        pids = new double[goodSeqs][goodSeqs];
        gaps = new double[goodSeqs][goodSeqs]; // gaps (and potentially mismatches)
        isNotUnique = new boolean[goodSeqs];
        // calculate id matrix and mpi
        for (int k = 25; k != OutAln[0].length() - 1; k++) { // -1 avoids line break char, 25 is clustal seq ID length
            lines:
            for (int i = 0; i != goodSeqs; i++) {
                double sum=0.0;
                // initiate gaps[] and pids[] to 0 ???????????????????????
                if (isNotUnique[i])
                    continue lines;
                for (int j = i + 1; j != goodSeqs; j++) {
                    try {

                        if (OutAln[i].charAt(k) == OutAln[j].charAt(k)) {
                            // this DP matrix makes shit easy!
                            if (OutAln[i].charAt(k) == 'A' || OutAln[i].charAt(k) == 'T'
                                    || OutAln[i].charAt(k) == 'C' || OutAln[i].charAt(k) == 'G'
                                    || OutAln[i].charAt(k) == 'U') { // U is just in case alignments are RNA
                                pids[i][j]++;
                                pids[j][i]++;
                            }
                        }
                        // this ignores "-:-"
                        else if (OutAln[i].charAt(k) != OutAln[j].charAt(k)) {
                            pids[j][i]++; // mismatch
                            if (OutAln[i].charAt(k) == '-' || OutAln[i].charAt(k) == '.')
                                gaps[i][j]++; // gap
                            if (OutAln[j].charAt(k) == '-' || OutAln[j].charAt(k) == '.')
                                gaps[j][i]++; // gap
                        }
                    } catch (Exception E) {
                        E.printStackTrace();
                        System.err.println("Caught Exception!\n");
                        System.err.println(i + " " + j + " " + k + " " + OutAln.length + " " + goodSeqs
                                + " FilteredTab[i]=" + FilteredTab[i].length() + " FilteredTab[j]=" + FilteredTab[j].length());
                        System.err.println("OutAln[i]=" + OutAln[i].length() + " OutAln[j]=" + OutAln[j].length()
                                + "\n" + OutAln[i] + "\n" + OutAln[j]);
                    }
                    // keep unique seqs ignoring gaps
                    if (k == OutAln[0].length() - 2) {
                        if (pids[j][i] - gaps[i][j] == pids[i][j] ||
                                pids[j][i] - gaps[j][j] == pids[i][j]) {
                            //both sequences are identical without gaps
                            //keep the longer one
                            if (gaps[i][j] > gaps[j][i])
                                isNotUnique[i] = true;
                            else
                                isNotUnique[j] = true; // this should also consider identical seqs
                        }

                            uniqueComps++;
                            // old mean pairwise identity ( considers gaps )
                            mpi = mpi + 100 * pids[i][j] / pids[j][i];

                            // classical average identity
                            mpi2 = mpi2 + 100 * pids[i][j] / Math.min(OutAln[i].replaceAll("[^ATCGU]", "").length(),
                                    OutAln[j].replaceAll("[^ATCGU]", "").length());


                    }
                }
            }

        }




        // calculate gaps, GC, shanon entropy, and Reverse Complement
        for (int k = 25; k != OutAln[0].length(); k++) {
            chars = new double[5];
            for (int i = 0; i != goodSeqs; i++) {
                /*if (isNotUnique[i]) {
                    if (k == OutAln[0].length() - 2)
                        uniqueSeqs--;
                    continue;
                }*/

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
            for (int z = 0; z != 5; z++)
                shanon = (chars[z] == 0) ? shanon + 0 : shanon + chars[z] / uniqueSeqs * (Math.log(chars[z] / uniqueSeqs) / Math.log(2));
        }
       // System.out.println( uniqueSeqs +"\t"+goodSeqs+"\t"+outCols+"\t"+totalChars[4]+"\t"+( outCols * goodSeqs));
        stats[0] = mpi / uniqueComps;                                                                        // Mean Pairwise ID
        stats[5] = mpi2 / uniqueComps;                                                                    // classical MPI
        for (int seq1 = 0; seq1 != goodSeqs; seq1++)
            for (int seq2 = seq1 + 1; seq2 != goodSeqs; seq2++)
                if (!isNotUnique[seq1] && !isNotUnique[seq2])
                    var = var + (double) Math.pow((100 * pids[seq1][seq2] / pids[seq2][seq1]) - stats[0], 2);
        stats[1] = var / uniqueComps;                                                                // Variance
        stats[2] = -1 * shanon / ((double) outCols);                                                        // Normalized Shanon entropy
        stats[3] = 100 * (totalChars[2] + totalChars[3]) / (totalChars[0] + totalChars[1] + totalChars[2] + totalChars[3]);       // GC content
        stats[4] = 100 * totalChars[4] / (outCols * goodSeqs);                                          // GAP content
/*

        //*********************************************************************
        //					R-scape scan & parse						*
        //*********************************************************************
        ///***************** 	R-scape scan & parse		******************
        // Write Sequences to Stockholm Format
        File dir = new File(dirPath + "/R-Scape");
        if (!dir.exists()){
            dir.mkdir();
        }
        String stkName =  key[1] + "_" + key[2] + ".stk";
        File stkFile = new File(dir + "/" + stkName);


        try {
            BufferedWriter WriteStockholm = new BufferedWriter(new FileWriter( stkFile ));
            WriteStockholm.write("# STOCKHOLM 1.0\n") ;

            for (int i =0; i < motifs.size(); i++){
                WriteStockholm.write(motifs.get(i)[0] + "\t"+motifs.get(i)[1] +"\n");
            }

            String gcRef = "#=GC RF "+ key[6] + "\n";

            String  gcSScon =  "#=GC SS_cons " +key[7] + "\n";

            WriteStockholm.write(gcRef);
            WriteStockholm.write(gcSScon);
            WriteStockholm.write("//" +"\n");


            WriteStockholm.close() ;
        } catch (IOException Fuck) {
            if (VERBOSE)
                System.err.println("couldn't write stockholm file");
            Fuck.printStackTrace();
            stkFile.delete() ;
            return;
        }


        try {
            rScape(stkFile, RSCAPEBINARY, dirPath);
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }



        String endFileName = ".fold.power";
        File f = new File(dirPath+ "/R-Scape");
        File[] matchingFiles = f.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return  name.endsWith(endFileName);
            }
        });

        String totalBasePairs = "";
        String expectedCovary = "";
        String observedCovary = "";

        try {
            Scanner scanner = new Scanner(matchingFiles[0]);

            //now read the file line by line...
            int lineNum = 0;
            int counter = 0;
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();
                lineNum++;

                if (line.startsWith("# BPAIRS expected to covary")) {
                    expectedCovary = line;

                } else if (line.startsWith("# BPAIRS observed to covary")) {
                    observedCovary = line;

                }else if (counter == 5) {
                    totalBasePairs = line;
                    counter ++;

                } else if (line.startsWith("#")) {
                    counter++;

                }
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        }
        String [] expCov= expectedCovary.split(" ");
        String [] totalBP = totalBasePairs.split(" ");
        String [] obsCov = observedCovary.split(" ");
        double covExp = Double.parseDouble(expCov[5]);
        double totalBasePair = Double.parseDouble(totalBP[totalBP.length - 1 ]);
        double observCov = Double.parseDouble(obsCov[obsCov.length - 1 ]);
        double percentAlignPower;
        double bpCovary;
        if (totalBasePair != 0.0) {
            percentAlignPower = covExp / totalBasePair * 100;
            bpCovary = observCov / totalBasePair * 100;
        } else {
            percentAlignPower =0.0;
            bpCovary =0.0;
        }

       //Delete file with R-scape info
        String[] entries = dir.list();
        for(String s:entries){
            File currentFile = new File (dir.getPath(),s);
            currentFile.delete();
        }
*/


        // save BED coords
        if (VERBOSE)
            System.out.println("- -> Calculating BED coords ");
        String BedFile = key[0] + "\t";
        BedFile = BedFile + key[1] + "\t" + key[2] + "\t";

        BedFile = BedFile + (int) goodSeqs + ":" + ((double) (int) (10 * stats[0]) / 10) + ":"      // MPI
                + ((double) (int) (10 * stats[5]) / 10) + ":"                      // CLASSIC MPI
                + ((double) (int) (10 * stats[4]) / 10) + ":"                     // GAPS
                + ((double) (int) (10 * Math.sqrt(stats[1])) / 10) + ":"            // STDEV
                + ((double) (int) (100 * stats[2]) / 100) + ":"                // SHANON
                + ((double) (int) (10 * stats[3]) / 10);                  //      GC
               // + (double) (int) percentAlignPower + ":"
               // + (double) (int) bpCovary + ":"
              //  + (double) (int) totalBasePair;
        if (VERBOSE)
            System.out.println("Pre SISSIz bed file: \n" + " " + BedFile);

        int random = (int) ((double) 10000 * Math.random());
        File Aln = new File(Path + "/" + BedFile.replaceAll("\t", "_") + ".aln." + random),    //
                AlnRC = new File(Path + "/" + BedFile.replaceAll("\t", "_") + "rc.aln." + random);  //
        // v v v v v v v v    INCLUSION STATS     v v v v v v v v v v v v v
        if (outCols > (FilteredTab[0].length()) / 2 && stats[4] <= 75 && stats[0] > 60) {
            //totalBasePair != 0.0 && !(percentAlignPower > 10 && bpCovary < 2)

            // Write Sequences to ALN Format
            try {
                BufferedWriter WriteClustal = new BufferedWriter(new FileWriter( Aln )),
                        WriteClustalRC = new BufferedWriter(new FileWriter( AlnRC ));
                WriteClustal.write("CLUSTAL format sucks\n\n") ;
                WriteClustalRC.write("CLUSTAL format sucks\n\n") ;
                for ( int y = 0 ; y != goodSeqs ; y++ ) {
                   // if ( !isNotUnique[ y ] ) {
                        WriteClustal.write( OutAln[ y ] ) ;
                        OutAlnRC[ y ] = OutAln[ y ].substring(0,25)+ OutAlnRC[ y ]  ;
                        WriteClustalRC.write( OutAlnRC[ y ] ) ;
                  //  }
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
            }
        } else {
            if (VERBOSE) {
                System.out.println("---> rejected alignment");
                System.out.println("     outcols = " + outCols + "\tuniqueseqs = " + uniqueSeqs +
                        "\tGAPS = " + stats[4] + "\n    PID = " + stats[0]);
                if (stats[0] < 5)
                    System.out.println("-----> SUPER LOW PID");
            }
            Aln.delete();
            AlnRC.delete();
            return;
        }
        String FinalBedFile = "",
                FinalBedFileRC = "",
                Antisense = (key[3].equals("+"))? "-" : "+";


//***************** 	SISSIz scan & parse		******************
        String[] SissizOutTab = new String[12];
        try {
            SissizOutTab = ScanSSZ(Path, BedFile, random);
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
        // delete empty alignments
        if (SissizOutTab == null || SissizOutTab[10] == null) {
            Aln.delete();
        } else {
            FinalBedFile = BedFile + ":" + SissizOutTab[1] + "_" + (int) (Double.parseDouble(SissizOutTab[10]) * -100) + "_" + key[3];
            // delete low scoring alignments
            if (Double.parseDouble(SissizOutTab[10]) > SSZR_THRESHOLD) {
                Aln.delete();
                if (PRINTALL) {
                    System.out.println(FinalBedFile.replaceAll("_", "\t"));
                }
            } else {
                //write bed and rename alignment
                System.out.println(FinalBedFile.replaceAll("_", "\t"));
                File NewFile = new File(Path + "/" + FinalBedFile.replaceAll("\t", "_") + ".aln");
                int file_count = 0;
                while (NewFile.exists()) {
                    file_count++;
                    NewFile = new File(Path + "/" + FinalBedFile.replaceAll("\t", "_") + ".aln_" + file_count);
                }
                boolean result = Aln.renameTo(NewFile);
            }
        }
        // * * * * * *  now for the RC  * * * * * *
        try {
            SissizOutTab = ScanSSZ(Path, BedFile + "rc", random);
            if (SissizOutTab == null) {
                AlnRC.delete();
            }
        } catch (IOException Fuck) {
            Fuck.printStackTrace();
            System.err.println("ScanSSZ failed in RC with ");
            for (int y = 0; y != goodSeqs; y++) {
                if (!isNotUnique[y]) {
                    System.err.println(OutAln[y]);
                }
            }
            AlnRC.delete();
        }
        if (SissizOutTab == null || SissizOutTab[10] == null) {
            AlnRC.delete();
        } else {
            FinalBedFileRC = BedFile + ":" + SissizOutTab[1] + "_" + (int) (Double.parseDouble(SissizOutTab[10]) * -100) + "_" + Antisense;
            // delete low scoring alignments
            if (Double.parseDouble(SissizOutTab[10]) > SSZR_THRESHOLD) {
                AlnRC.delete();
                if (PRINTALL) {
                    System.out.println(FinalBedFileRC.replaceAll("_", "\t"));
                }
            } else {
                //write bedRC and rename alignment
                System.out.println(FinalBedFileRC.replaceAll("_", "\t"));
                File NewFile = new File(Path + "/" + FinalBedFileRC.replaceAll("\t", "_") + ".aln");
                int file_count = 0;
                while (NewFile.exists()) {
                    file_count++;
                    NewFile = new File(Path + "/" + FinalBedFileRC.replaceAll("\t", "_") + ".aln_" + file_count);
                }
                boolean result = AlnRC.renameTo(NewFile);
            }
            return;
        }
       // Aln.delete() ;
       // AlnRC.delete() ;


    }


//*********************************************************************
    //					Rscape scan & parse						*
    //*********************************************************************


    private void rScape(File stockholmFile, String RSCAPEBINARY, String path)throws IOException,
            InterruptedException {
        String pathway = path + "/R-Scape";

        String cmd = RSCAPEBINARY + " --fold" + " --nofigures "
                + stockholmFile;
        Process process = Runtime.getRuntime().exec(cmd, null, new File(pathway));
        BufferedReader reader =
                new BufferedReader(new InputStreamReader(process.getInputStream()));
        while ((reader.readLine()) != null) {
        }
        process.waitFor();
    }
        //*********************************************************************
        //					SISSIz scan & parse						*
        //*********************************************************************
        // sissiz-di       cluster.109999_step.aln  8       150     0.8759  0.8542  0.0094  -13.88  -8.20   3.48    -1.63
        protected static String[] ScanSSZ (String Path, String BedFile, int id ) throws
        IOException {
            //stats[0] Mean Pairwise ID
            //stats[1] Variance
            //stats[2] Normalized Shanon entropy
            //stats[3] GC content
            //stats[4] GAP content
            //stats[5] classical MPI
            String SissizOutTab[] = new String[12];
            String Output = "", Error = "";
            String Command = SSZBINARY;

            Command = Command + " -j " + Path + "/" + BedFile.replaceAll("\t", "_") + ".aln." + id; // RIBOSUM scoring

            try {
                long now = System.currentTimeMillis();
                long timeoutInMillis = 1000L * 300;                          // max 5 minutes
                long finish = now + timeoutInMillis;
                // launch initial SISSIz call
                Process Sissiz = Runtime.getRuntime().exec(Command);
                BufferedReader SissizErr = new BufferedReader(new InputStreamReader(Sissiz.getErrorStream()));
                if (VERBOSE)
                    System.out.println(": Running " + Command);
                while (isAlive(Sissiz)) {
                    Thread.sleep(100);
                    if (System.currentTimeMillis() > finish) {
                        if (VERBOSE)
                            System.out.println("SISSIz failed to run within time :-(");
                        SissizErr.close();
                        Sissiz.destroy();
                        return null;


                        }
                }
                SissizErr.close();
                // get Output if process didn't complete in recursion
                if (SissizOutTab[0] == null) {
                    BufferedReader SissizOut = new BufferedReader(new InputStreamReader(Sissiz.getInputStream()));
                    while ((Output = SissizOut.readLine()) != null) {
                        if (Output.length() > 6 && Output.substring(0, 6).equals("sissiz")) {
                            if (VERBOSE)
                                System.out.println(Output);
                            SissizOutTab = Output.split("\\s");
                            SissizOutTab[1] = "r";
                        }
                    }
                    SissizOut.close();
                }
                // rerun SISSIz if output is dodgy
                try {
                    /*>>>>>>*/
                    if (Double.parseDouble(SissizOutTab[7]) == 0
                            || (Math.abs(Double.parseDouble(SissizOutTab[9])) < 0.5)){ // variance of simulated alignments
                        if (VERBOSE)
                            System.out.println("SISSIz gave dodgy output... retrying with RIBOSUM");
                        SissizOutTab = new String[12];
                        now = System.currentTimeMillis();
                        finish = now + timeoutInMillis;
                        Sissiz = Runtime.getRuntime().exec(SSZBINARY + " -j " + Path + "/" + BedFile.replaceAll("\t", "_") + ".aln." + id);
                        if (VERBOSE)
                            System.out.println("Running " + SSZBINARY + " -j " + Path + "/" + BedFile.replaceAll("\t", "_") + ".aln." + id);
                        SissizErr = new BufferedReader(new InputStreamReader(Sissiz.getErrorStream()));
                        while (isAlive(Sissiz)) {
                            Thread.sleep(100);  // do we need this?
                            if (System.currentTimeMillis() > finish) {
                                if (VERBOSE)
                                    System.out.println("SISSIz failed to run within time");
                                SissizErr.close();
                                Sissiz.destroy();
                                return null;
                            } else {
                                // avoid infinity loop
                                if (SissizErr.ready()) {
                                    Error = SissizErr.readLine();
/*								if ( Error.length() > 17 && Error.substring(0,17).equals( "WARNING: Negative") && counter <2 ) {
									if (VERBOSE) {
										System.out.println( Error ) ;
										System.out.println("         Adding extra flanking sites");
									}
									SissizErr.close();
									Sissiz.destroy();
									// launching with more flanking sites
									SissizOutTab = ScanSSZ( Path, BedFile, stats, 4, id );
									break;
								}
*/
                                }
                            }
                        }
                        SissizErr.close();
                        if (SissizOutTab[0] == null) {
                            BufferedReader SissizOut = new BufferedReader(new InputStreamReader(Sissiz.getInputStream()));
                            while ((Output = SissizOut.readLine()) != null) {
                                if (Output.length() > 6 && Output.substring(0, 6).equals("sissiz")) {
                                    if (VERBOSE)
                                        System.out.println(Output + "  ...DONE !!\n");
                                    SissizOutTab = Output.split("\\s");
                                    SissizOutTab[1] = "r";
                                }
                            }
                            SissizOut.close();
                        }
                    }

                } catch (Exception Fuck) {
                    System.err.println("  with file = " + BedFile.replaceAll("\t", "_") + ".aln");
                    if (SissizOutTab[7] == null)
                        System.err.println("Null output data");
                    Fuck.printStackTrace();
                }
            } catch (Exception err) {
                System.out.println(" Caught Error!\n ----> " + Command + "\n  counter--> " );
                System.err.println("!!!Caught Error!\n ----> " + Command + "\n  counter--> " );
                //String [] OutAln = new String[ goodSeqs ] ;
                err.printStackTrace();
                System.err.println("===============================");
            }
            return SissizOutTab;
        }

    //*********************************************************************
    //						Sample process						*
    //*********************************************************************
    private static boolean isAlive( Process p ) {
        try {
            p.exitValue();
            return false;
        } catch (IllegalThreadStateException e) {
            return true;
        }
    }
    public void setSszR(double newValue){
        SSZR_THRESHOLD = newValue;
    }
}
