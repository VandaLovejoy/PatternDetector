import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class AssociativeList {
    private ArrayList<String []>associativeList;
    int numbZeros;
    // Constructor
    public AssociativeList() { associativeList = new ArrayList<>();}

    public HashMap<String[], ArrayList> addSpecies(String pathFile, int blockAln) throws IOException {
        HashMap<String[], ArrayList> motifs = new HashMap<>();
        String lastDigit = "";
        String finalName = "";
        lastDigit += String.valueOf(blockAln);
        numbZeros = 4 - lastDigit.length();
        int i = 0;
        while (i < numbZeros) {
            finalName += String.valueOf(0);
            i++;
        }
        finalName += lastDigit;

        File file = new File(pathFile + "alifold_" + finalName + ".stk");

        BufferedReader reader = new BufferedReader(new FileReader(file));
        String currentLine = reader.readLine();
        readBlocks:
        while (currentLine != null ) {

            if (currentLine.length() != 2 && currentLine.substring(0, 7).equals("#=GF ID")) {
                String[] coordinates = currentLine.split("_");
                int begCord = Integer.parseInt(coordinates[4]);
                int endCord = Integer.parseInt(coordinates[3]);
                int lengthSeq = begCord - endCord;
                if (lengthSeq >= 50) {
                    associativeList = new ArrayList<>();
                    String[] coordinates2 = new String[2];
                    coordinates2[0]= coordinates[3];
                    coordinates2[1]= coordinates[4];
                    currentLine = reader.readLine();
                    currentLine = reader.readLine();

                    while (currentLine.charAt(0) != '#') {
                        String[] species = currentLine.split(" ", 2);
                        species[1]=species[1].trim();

                        associativeList.add(species);

                        currentLine = reader.readLine();
                    }
                    motifs.put(coordinates2, associativeList);
                    String[] endMotif = new String[]{"//"};
                   //  associativeList.add(endMotif);
                     currentLine =reader.readLine();
                     continue readBlocks;
                } else {
                   currentLine =reader.readLine();
                    continue readBlocks;
                }
            } else {
                currentLine =reader.readLine();
                continue readBlocks;
            }

        }
        reader.close();
        file.delete();
        return motifs;
    }

    public char[] getAlnTab(String [] seq){
        char[] AlnTab = new char[seq[0].length()];
        if (seq.length>1 && seq[1].length()> 10) {
            AlnTab = seq[1].toCharArray();
        }
        return AlnTab;

    }

}
