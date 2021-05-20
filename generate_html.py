import lrst_global


class ST_HTML_Generator:
   def __init__(self, para_list):
      self.image_key_list = para_list[0]
      self.header_info = para_list[1];
      self.input_para = para_list[2];
      if len( self.input_para["input_files"] ) > 1:
         self.more_input_files = True;
      else: self.more_input_files = False;      

   def generate_header( self):
      self.html_writer = open( self.input_para["output_folder"]+'/'+self.input_para["out_prefix"]+"statistics.html", 'w')
      self.html_writer.write("<html>")
      self.html_writer.write( "<head>" )
      self.html_writer.write("<title>" )
      self.html_writer.write( self.header_info )
      self.html_writer.write("</title>" )

      self.html_writer.write('''<style  type="text/css">
 @media screen {
  div.summary {
    width: 18em;
    position:fixed;
    top: 3em;
    margin:1em 0 0 1em;
  }
  
  div.main {
    display:block;
    position:absolute;
    overflow:auto;
    height:auto;
    width:auto;
    top:4.5em;
    bottom:2.3em;
    left:18em;
    right:0;
    border-left: 1px solid #CCC;
    padding:0 0 0 1em;
    background-color: white;
    z-index:1;
  }
  
  div.header {
    background-color: #EEE;
    border:0;
    margin:0;
    padding: 0.5em;
    font-size: 200%;
    font-weight: bold;
    position:fixed;
    width:100%;
    top:0;
    left:0;
    z-index:2;
  }

  div.footer {
    background-color: #EEE;
    border:0;
    margin:0;
	padding:0.5em;
    height: 1.3em;
	overflow:hidden;
    font-size: 100%;
    font-weight: bold;
    position:fixed;
    bottom:0;
    width:100%;
    z-index:2;
  }
  
  img.indented {
    margin-left: 3em;
  }
 }
 
 @media print {
	img {
		max-width:100% !important;
		page-break-inside: avoid;
	}
	h2, h3 {
		page-break-after: avoid;
	}
	div.header {
      background-color: #FFF;
    }
	
 }
 
 body {    
  font-family: sans-serif;   
  color: #000;   
  background-color: #FFF;
  border: 0;
  margin: 0;
  padding: 0;
  }
  
  div.header {
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 200%;
  font-weight: bold;
  width:100%;
  }    
  
  #header_title {
  display:inline-block;
  float:left;
  clear:left;
  }
  #header_filename {
  display:inline-block;
  float:right;
  clear:right;
  font-size: 50%;
  margin-right:2em;
  text-align: right;
  }

  div.header h3 {
  font-size: 50%;
  margin-bottom: 0;
  }
  
  div.summary ul {
  padding-left:0;
  list-style-type:none;
  }
  
  div.summary ul li img {
  margin-bottom:-0.5em;
  margin-top:0.5em;
  }
	  
  div.main {
  background-color: white;
  }
      
  div.module {
  padding-bottom:1.5em;
  padding-top:1.5em;
  }
	  
  div.footer {
  background-color: #EEE;
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 100%;
  font-weight: bold;
  width:100%;
  }


  a {
  color: #000080;
  }

  a:hover {
  color: #800000;
  }
      
  h2 {
  color: #800000;
  padding-bottom: 0;
  margin-bottom: 0;
  clear:left;
  }

  table { 
  margin-left: 3em;
  text-align: center;
  }
  
  th { 
  text-align: center;
  background-color: #000080;
  color: #FFF;
  padding: 0.4em;
  }      
  
  td { 
  font-family: monospace; 
  text-align: left;
  background-color: #EEE;
  color: #000;
  padding: 0.4em;
  }

  img {
  padding-top: 0;
  margin-top: 0;
  border-top: 0;
  }

  
  p {
  padding-top: 0;
  margin-top: 0;
  }
      </style>''')

      self.html_writer.write( "</head>" )
      self.html_writer.write( "<body>" )  
      self.html_writer.write( '''<div class="header"> 
         <div id="header_title">LongReadDS Report</div>
         <div id="header_filename">
             <script> document.write(new Date().toLocaleDateString()); </script>
      ''')
      # for _af in self.input_para["input_files"]:
      #   self.html_writer.write( "<br/>"+_af); 
      self.html_writer.write( "<br/>"+ self.input_para["input_files"][0] )
      self.html_writer.write( '''       
         </div>
      </div>''')
 
   def generate_left( self ):
      self.html_writer.write('<div class="summary">');
      self.html_writer.write('<h2>Summary</h2>')
      self.html_writer.write('<ul>')
 
      _imki = 0;     
      for _imk in self.image_key_list:
          self.html_writer.write('<li>')
          self.html_writer.write('<a href="#lrst'+str(_imki)+'">'+lrst_global.plot_filenames[_imk]['title']+'</a>')
          _imki += 1;
          self.html_writer.write('</li>')
      if self.more_input_files:
          self.html_writer.write('<li>')
          self.html_writer.write('<a href="#lrst'+str(_imki)+'">Input files</a>')
          _imki += 1;
          self.html_writer.write('</li>')

      self.html_writer.write("</ul>")
      self.html_writer.write('</div>')

   def generate_right( self ):
      self.html_writer.write('<div class="main">')
      _imki = 0;
      for _imk in self.image_key_list:
         self.html_writer.write('<div class="module">');
         self.html_writer.write('<h2 id="lrst'+str(_imki)+'">'+lrst_global.plot_filenames[_imk]['description']+'</h2><p>')
         self.html_writer.write('<img class="indented" src="'+lrst_global.plot_filenames[_imk]['file']+'" alt="'+lrst_global.plot_filenames[_imk]['description']+'" width="800" height="600"/></p>')
         self.html_writer.write('</div>')
         _imki += 1;
      if self.more_input_files:
         self.html_writer.write('<div class="module">');
         self.html_writer.write('<h2 id="lrst'+str(_imki)+'">The list of input files: '+str(len( self.input_para["input_files"] ))+'</h2><p>')
         for _af in self.input_para["input_files"]:
            self.html_writer.write( "<br/>"+_af);
         self.html_writer.write('</p></div>')
         _imki += 1;

      self.html_writer.write('</div>')

   def generate_end( self ):
      self.html_writer.write( '<div class="footer"> Generated by <a href="https://github.com/WGLab/LongReadDS">LongReadDS</a> </div>' )

      self.html_writer.write( "</body>" )
      self.html_writer.write("</html>")
      self.html_writer.close();

   def generate_st_html( self ):
      self.generate_header( );
   
      self.generate_left( );
      self.generate_right( );

      self.generate_end( );


