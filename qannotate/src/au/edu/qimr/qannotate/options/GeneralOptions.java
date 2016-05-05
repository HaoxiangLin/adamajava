package au.edu.qimr.qannotate.options;

import static java.util.Arrays.asList;

import java.awt.List;

import au.edu.qimr.qannotate.Messages;
import au.edu.qimr.qannotate.options.Options.MODE;
import joptsimple.OptionSet;

/*
 * parse command line to options. 
 */
public class GeneralOptions extends Options {
    private int bufferSize = 0;


    /**
     * check command line and store arguments and option information
     */
 	
 	// public Options(){  parser = new OptionParser(); this.Mode = null; } 
 	public GeneralOptions(Options.MODE mode){ 		super(mode);   }
 	
    @Override
    public boolean parseArgs(final String[] args) throws Exception{  	
    	 
        parser.acceptsAll( asList("h", "help"), Messages.getMessage("HELP_OPTION_DESCRIPTION"));
        parser.acceptsAll( asList("i", "input"), Messages.getMessage("INPUT_DESCRIPTION")).withRequiredArg().ofType(String.class).describedAs("input vcf");
        parser.acceptsAll( asList("o", "output"), Messages.getMessage("OUTPUT_DESCRIPTION")).withRequiredArg().ofType(String.class).describedAs("output vcf"); 
        parser.acceptsAll( asList("d", "database"), Messages.getMessage("DATABASE_DESCRIPTION")).withRequiredArg().ofType(String.class).describedAs("database file"); 
        parser.accepts("mode", "run mode").withRequiredArg().ofType(String.class).describedAs("mode");
        parser.accepts("log", LOG_DESCRIPTION).withRequiredArg().ofType(String.class);
        parser.accepts("loglevel",  LOG_LEVEL_OPTION_DESCRIPTION).withRequiredArg().ofType(String.class);
        if(Mode.equals(MODE.trf))
            parser.accepts("buffer", "check TRF region on both side of indel within this nominated size" ).withRequiredArg().ofType(Integer.class);//.describedAs("integer");
      
        OptionSet options = parser.parse(args);            
         
        if(options.has("h") || options.has("help")){        	
        	displayHelp(Mode);
            return false;
        }
               
        if( !options.has("log")){
            System.out.println(Messages.getMessage("LOG_OPTION_DESCRIPTION"));            
            return false;
        } else{  
        	logFileName = (String) options.valueOf("log");  	
        	logLevel = (String) options.valueOf("loglevel");
        }   
        
        commandLine = Messages.reconstructCommandLine(args) ;
        
        //check IO
        inputFileName = (String) options.valueOf("i") ;      	 
        outputFileName = (String) options.valueOf("o") ; 
        databaseFileName = (String) options.valueOf("d") ; 
        if(options.has("buffer"))
        	bufferSize = (Integer) options.valueOf("buffer");   
        			
        String[] inputs = new String[]{ inputFileName,databaseFileName} ;
        String[] outputs = new String[]{outputFileName};
        String [] ios = new String[inputs.length + outputs.length];
        System.arraycopy(inputs, 0, ios, 0, inputs.length);
        System.arraycopy(outputs, 0, ios, inputs.length, outputs.length);
        return checkInputs(inputs )  && checkOutputs(outputs ) && checkUnique(ios);
       
    } 
    
    public int getBufferSize(){
    	
    	 if(Mode.equals(MODE.trf))
    		 return bufferSize;
    	 
    	 return -1; 
    }
   
}