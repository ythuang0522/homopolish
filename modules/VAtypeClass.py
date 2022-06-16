class FixSNP():
    def __init__(self):
        self._seq = ""
        self._reads_file = ""
        self._sibilings = ""
        self._sketch_path = ""
        self._output_dir = ""
        self._mash_threshold = 0.95
        self.dl_contig_nums = 20
        self._genome_file = ""
        self._contig_id = ""
        self._rmFilePath = "Homo"   
        self._fix_genome_file = ""
        self._true_genome_file = ""
        self._getFixPosCSV = True
        self._getMissPosCSV = False
        self._getErrorPosCSV = False
        self._spPattern = ""
        self._bam_file = ""
        self._thread = 16
        self._gen_description = ""
       
       


    @property
    def gen_desc(self):
        return self._gen_description

    @gen_desc.setter
    def gen_desc(self, value):
        self._gen_description = value
               
        
    @property
    def thread(self):
        return self._thread

    @thread.setter
    def thread(self, value):
        self._thread = value
        
    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, value):
        self._seq = value
        
    @property
    def bamFile(self):
        return self._bam_file

    @bamFile.setter
    def bamFile(self, value):
        self._bam_file = value
        
    

    @property
    def spPattern(self):
        return self._spPattern

    @spPattern.setter
    def spPattern(self, value):
        self._spPattern = value

    
    #reads file
    @property
    def reads_file(self):
        return self._reads_file

    @reads_file.setter
    def reads_file(self, value):
        self._reads_file = value
        
        
        
#   sibilings:draft_gem的 sibilings          
    @property
    def sib_files(self):
        return self._sibilings

    @sib_files.setter
    def sib_files(self, value):
        self._sibilings = value
       
        
    #sketch path         
    @property
    def sketch_path(self):
        return self._sketch_path

    @sketch_path.setter
    def sketch_path(self, value):
        self._sketch_path = value
        
    #sketch path         
    @property
    def output_dir(self):
        return self._output_dir

    @output_dir.setter
    def output_dir(self, value):
        self._output_dir = value

    #mash_threshold  相似度     
    @property
    def mash_threshold(self):
        return self._mash_threshold

    @mash_threshold.setter
    def mash_threshold(self, value):
        self._mash_threshold = value
        
    #設定下載siblings數量    
    @property
    def dl_contig_nums(self):
        return self._dl_contig_nums

    @dl_contig_nums.setter
    def dl_contig_nums(self, value):
        self._dl_contig_nums = value
        
    #給mash的fa檔案   
    @property
    def draft_genome_file(self):
        return self._genome_file

    @draft_genome_file.setter
    def draft_genome_file(self, value):
        self._genome_file = value
        
    #fa檔案的contig_id  
    @property
    def contig_id(self):
        return self._contig_id

    @contig_id.setter
    def contig_id(self, value):
        self._contig_id = value
        
    
    #remove siblings gz file path  
    @property
    def rmFilePath(self):
        return self._rmFilePath

    @rmFilePath.setter
    def rmFilePath(self, value):
        self._rmFilePath = value
          
    #fix genome file
    @property
    def fix_genome_file(self):
        return self._fix_genome_file

    @fix_genome_file.setter
    def fix_genome_file(self, value):
        self._fix_genome_file = value
        
        
    #true genome file
    @property
    def true_genome_file(self):
        return self._true_genome_file

    @true_genome_file.setter
    def true_genome_file(self, value):
        self._true_genome_file = value   
                
    
    #flag of fixPosCSV
    @property
    def get_fixCSV_Flag(self):
        return self._getFixPosCSV

    @get_fixCSV_Flag.setter
    def get_fixCSV_Flag(self, value):
        self._getFixPosCSV = value   
        
        
    #flag of MissPosCSV
    @property
    def get_MissCSV_Flag(self):
        return self._getMissPosCSV

    @get_MissCSV_Flag.setter
    def get_MissCSV_Flag(self, value):
        self._getMissPosCSV = value   
    
    #flag of ErrorPosCSV
    @property
    def get_EorCSV_Flag(self):
        return self._getErrorPosCSV

    @get_EorCSV_Flag.setter
    def get_EorCSV_Flag(self, value):
        self._getErrorPosCSV = value   
        
        
  
    
    
        
