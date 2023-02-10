"""
generate_html.py: Generate the HTML file from our plot images.
"""

import base64
from src import lrst_global  # Contains our image filepaths


class ST_HTML_Generator:
    def __init__(self, para_list, static=True):
        self.image_key_list = para_list[0]  # List of statistics variables to look for
        self.header_info = para_list[1]  # The webpage's title
        self.input_para = para_list[2]  # List of the input parameters used for these statistics
        self.static = static  # Static vs. dynamic webpage boolean

        if len(self.input_para["input_files"]) > 1:
            self.more_input_files = True;
        else:
            self.more_input_files = False;

    def generate_header(self):
        if self.static:
            self.html_writer = open(
                self.input_para["output_folder"] + '/' + self.input_para["out_prefix"] + "statistics.html", 'w')
        else:
            self.html_writer = open(
                self.input_para["output_folder"] + '/' + self.input_para["out_prefix"] + "statistics_dynamic.html", 'w')
        self.html_writer.write("<html>")
        self.html_writer.write("<head>")
        self.html_writer.write("<title>")
        self.html_writer.write(self.header_info)
        self.html_writer.write("</title>")
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
  
  li {
  margin: 10px 0;
  }
      </style>''')

        self.html_writer.write("</head>")
        self.html_writer.write("<body>")
        self.html_writer.write('''<div class="header"> 
         <div id="header_title">{} Report</div>
         <div id="header_filename">
             <script> document.write(new Date().toLocaleDateString()); </script>
      '''.format(lrst_global.prg_name))
        # for _af in self.input_para["input_files"]:
        #   self.html_writer.write( "<br/>"+_af);
        # self.html_writer.write( "<br/>"+ self.input_para["input_files"][0] )
        self.html_writer.write('''       
         </div>
      </div>''')

    def generate_left(self):
        self.html_writer.write('<div class="summary">');
        self.html_writer.write('<h2>Summary</h2>')
        self.html_writer.write('<ul>')

        _imki = 0
        for _imk in self.image_key_list:
            self.html_writer.write('<li>')
            self.html_writer.write(
                '<a href="#lrst' + str(_imki) + '">' + lrst_global.plot_filenames[_imk]['title'] + '</a>')
            _imki += 1;
            self.html_writer.write('</li>')

        self.html_writer.write('<li>')
        self.html_writer.write('<a href="#lrst' + str(_imki) + '">List of input files</a>')
        _imki += 1;
        self.html_writer.write('</li>')

        self.html_writer.write("</ul>")
        self.html_writer.write('</div>')

    def generate_right(self):
        self.html_writer.write('<div class="main">')
        _imki = 0
        for _imk in self.image_key_list:
            self.html_writer.write('<div class="module">')
            self.html_writer.write(
                '<h2 id="lrst' + str(_imki) + '">' + lrst_global.plot_filenames[_imk]['description'] + '</h2><p>')
            # self.html_writer.write('<img class="indented" src="'+lrst_global.plot_filenames[_imk]['file']+'"
            # alt="'+lrst_global.plot_filenames[_imk]['description']+'" width="600" height="450"/></p>')

            if 'dynamic' in lrst_global.plot_filenames[_imk] and self.static == False:
                self.html_writer.write(lrst_global.plot_filenames[_imk]['dynamic'])

            else:
                if _imk == "basic_st":
                    self.html_writer.write(lrst_global.plot_filenames["basic_st"]['detail'])
                else:
                    m_image_file = open(
                        self.input_para["output_folder"] + '/' + lrst_global.plot_filenames[_imk]['file'], 'rb');
                    self.html_writer.write('<img class="indented" src="data:image/png;base64,' + base64.b64encode(
                        m_image_file.read()).decode('utf-8') + '" alt="' + lrst_global.plot_filenames[_imk][
                                               'description'] + '" width="800" height="600"/></p>')
                    m_image_file.close()

            self.html_writer.write('</div>')

            _imki += 1

        self.html_writer.write('<div class="module">')
        self.html_writer.write('<h2 id="lrst' + str(_imki) + '">File count = ' + str(
            len(self.input_para["input_files"])) + '</h2><p>')
        for _af in self.input_para["input_files"]:
            self.html_writer.write("<br/>" + _af)
        self.html_writer.write('</p></div>')
        _imki += 1

        self.html_writer.write('</div>')

    def generate_left_signal_data(self, read_names):
        """Generate links in the left panel."""
        self.html_writer.write('<div class="summary">');
        self.html_writer.write('<h2>Summary</h2>')
        self.html_writer.write('<ul>')


        # Add the summary table section link
        url_index = 0
        self.html_writer.write('<li>')
        self.html_writer.write(
            '<a href="#lrst' + str(url_index) + '">' + lrst_global.plot_filenames["basic_st"]['title'] + '</a>')
        url_index += 1

        # Add the read plot section list
        for read_name in read_names:
            self.html_writer.write('<li>')
            self.html_writer.write(
                '<a href="#lrst' + str(url_index) + '">' + read_name + '</a>')
            url_index += 1

        # Add the input files section link
        self.html_writer.write('<li>')
        self.html_writer.write('<a href="#lrst' + str(url_index) + '">Input files</a>')
        url_index += 1
        self.html_writer.write('</li>')
        self.html_writer.write("</ul>")
        self.html_writer.write('</div>')

    def generate_right_signal_data(self, read_names, signal_plot):
        """Generate tables and plots in the right section."""

        self.html_writer.write('<div class="main">')

        # Add the summary table section
        url_index = 0
        self.html_writer.write('<div class="module">')
        self.html_writer.write(
            '<h2 id="lrst' + str(url_index) + '">' + lrst_global.plot_filenames["basic_st"]['description'] + '</h2><p>')
        self.html_writer.write(lrst_global.plot_filenames["basic_st"]['detail'])
        url_index += 1

        # Add the read plot section
        for read_name in read_names:
            self.html_writer.write('<div class="module">')

            # Set the description
            description_text = "Basecall signal"
            self.html_writer.write(
                '<h2 id="lrst' + str(url_index) + '">' + description_text + '</h2><p>')

            # Add the plot
            html_plot = signal_plot.get(read_name)
            self.html_writer.write(html_plot)
            self.html_writer.write('</div>')
            url_index += 1

        # Add the input files section
        self.html_writer.write('<div class="module">')
        self.html_writer.write('<h2 id="lrst' + str(url_index) + '">File count = ' + str(
            len(self.input_para["input_files"])) + '</h2><p>')
        for input_file_str in self.input_para["input_files"]:
            self.html_writer.write("<br/>" + input_file_str)
        self.html_writer.write('</p></div>')
        self.html_writer.write('</div>')

    def generate_end(self):
        self.html_writer.write(
            '<div class="footer"> Generated by <a href="https://github.com/WGLab/{}">{}</a> </div>'.format(
                lrst_global.prg_name, lrst_global.prg_name))

        self.html_writer.write("</body>")
        self.html_writer.write("</html>")
        self.html_writer.close()

    def generate_st_html(self, signal_plots=None):
        """
        Top-level function for generating the HTML.
        """
        if signal_plots is None:
            # Format base QC
            self.generate_header()
            self.generate_left()
            self.generate_right()
            self.generate_end()
        else:
            # Format signal QC
            self.generate_header()
            read_names = signal_plots.keys()
            self.generate_left_signal_data(read_names)
            self.generate_right_signal_data(read_names, signal_plots)
            self.generate_end()
