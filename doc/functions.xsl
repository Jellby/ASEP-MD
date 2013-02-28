<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:param name="lang">en</xsl:param>

<xsl:variable name="tag_default">
  <xsl:choose>
    <xsl:when test="$lang='en'"><xsl:text>Default</xsl:text></xsl:when>
    <xsl:when test="$lang='es'"><xsl:text>Predeterminado</xsl:text></xsl:when>
  </xsl:choose>
</xsl:variable>

<xsl:variable name="tag_allowed">
  <xsl:choose>
    <xsl:when test="$lang='en'"><xsl:text>Allowed values</xsl:text></xsl:when>
    <xsl:when test="$lang='es'"><xsl:text>Valores permitidos</xsl:text></xsl:when>
  </xsl:choose>
</xsl:variable>

<xsl:variable name="tag_other">
  <xsl:choose>
    <xsl:when test="$lang='en'"><xsl:text>any other value</xsl:text></xsl:when>
    <xsl:when test="$lang='es'"><xsl:text>cualquier otro valor</xsl:text></xsl:when>
  </xsl:choose>
</xsl:variable>

<xsl:template name="wrap-string" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:param name="str" />
  <xsl:param name="wrap-col" />
  <xsl:param name="break-mark" />
  <xsl:param name="indent" />
  <xsl:param name="first-indent" />
  <xsl:param name="pos" select="0" />
  <xsl:param name="start" select="1" />
  <xsl:choose>

    <!-- start of first line: write $first-indent -->
    <xsl:when test="$start">
      <xsl:copy-of select="$first-indent"/>
      <xsl:call-template name="wrap-string">
        <xsl:with-param name="str" select="$str" />
        <xsl:with-param name="wrap-col" select="$wrap-col" />
        <xsl:with-param name="break-mark" select="$break-mark" />
        <xsl:with-param name="indent" select="$indent" />
        <xsl:with-param name="pos" select="0" />
        <xsl:with-param name="start" select="0" />
      </xsl:call-template>
    </xsl:when>

    <!-- not the last word -->
    <xsl:when test="contains( $str, ' ' )">
      <xsl:variable name="first-word" select="substring-before( $str, ' ' )" />
      <xsl:variable name="pos-now" select="$pos + ( $pos > 0 ) + string-length( $first-word )" />
      <xsl:choose>

        <!-- a word that is not the first in a line and does not fit: write $break-mark and $indent -->
        <xsl:when test="$pos > 0 and $pos-now > $wrap-col">
          <xsl:copy-of select="$break-mark" />
          <xsl:copy-of select="$indent" />
          <xsl:call-template name="wrap-string">
            <xsl:with-param name="str" select="$str" />
            <xsl:with-param name="wrap-col" select="$wrap-col" />
            <xsl:with-param name="break-mark" select="$break-mark" />
            <xsl:with-param name="indent" select="$indent" />
            <xsl:with-param name="pos" select="0" />
            <xsl:with-param name="start" select="0" />
          </xsl:call-template>
        </xsl:when>

        <!-- other words: write a space if it's not the first of the line, then the word -->
        <xsl:otherwise>
          <xsl:if test="$pos > 0"><xsl:text> </xsl:text></xsl:if>
          <xsl:value-of select="$first-word" />
          <xsl:call-template name="wrap-string">
            <xsl:with-param name="str" select="substring-after( $str, ' ' )" />
            <xsl:with-param name="wrap-col" select="$wrap-col" />
            <xsl:with-param name="break-mark" select="$break-mark" />
            <xsl:with-param name="indent" select="$indent" />
            <xsl:with-param name="pos" select="$pos-now + ( $pos = 0 ) * ( $start = 0 ) * string-length($indent)" />
            <xsl:with-param name="start" select="0" />
          </xsl:call-template>
        </xsl:otherwise>

      </xsl:choose>
    </xsl:when>

    <!-- last word: write it in this line or in the next -->
    <xsl:otherwise>
      <xsl:choose>
        <xsl:when test="$pos + ( $pos > 0 ) + string-length( $str ) > $wrap-col">
          <xsl:copy-of select="$break-mark" />
          <xsl:copy-of select="$indent" />
          <xsl:value-of select="$str" />
        </xsl:when>
        <xsl:otherwise>
          <xsl:if test="$pos > 0"><xsl:text> </xsl:text></xsl:if>
          <xsl:value-of select="$str" />
        </xsl:otherwise>
      </xsl:choose>
    </xsl:otherwise>

  </xsl:choose>
</xsl:template>

<xsl:template name="write-block">
  <xsl:param name="text"/>
  <xsl:param name="first-indent"/>
  <xsl:param name="indent"/>
  <xsl:choose>
    <xsl:when test="contains( $text, '¬' )">
      <!-- a block between "¬{" and "¬}" is rendered as an equation-->
      <xsl:call-template name="write-block">
        <xsl:with-param name="text"><xsl:copy-of select="substring-before( $text, '¬{' )"/></xsl:with-param>
      </xsl:call-template>
      <xsl:call-template name="write-equation">
        <xsl:with-param name="text"><xsl:copy-of select="substring-before( substring-after( $text, '¬{' ), '¬}' )"/></xsl:with-param>
        <xsl:with-param name="wrap-col" select="80"/>
      </xsl:call-template>
      <xsl:call-template name="write-block">
        <xsl:with-param name="text"><xsl:copy-of select="substring-after( $text, '¬}' )"/></xsl:with-param>
        <xsl:with-param name="first-indent" select="$indent"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="wrap-string">
        <xsl:with-param name="str"><xsl:copy-of select="normalize-space($text)"/></xsl:with-param>
        <xsl:with-param name="wrap-col" select="80"/>
        <xsl:with-param name="first-indent" select="$first-indent"/>
        <xsl:with-param name="break-mark" select="$newline"/>
        <xsl:with-param name="indent" select="$indent"/>
      </xsl:call-template>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template name="write-equation">
  <xsl:param name="text"/>
  <xsl:param name="wrap-col" />
  <xsl:variable name="equation">
    <xsl:choose>
      <xsl:when test="contains( $text, '¬' )">
        <xsl:value-of select="normalize-space(substring-before( $text, '¬' ))"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of name="equation" select="normalize-space( $text )"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:variable>
  <xsl:variable name="label" select="normalize-space(substring-after( $text, '¬' ))"/>
  <xsl:variable name="eqlength" select="string-length( $equation )"/>
  <xsl:variable name="lalength" select="string-length( $label )"/>
  <xsl:variable name="spleft" select="floor ( ( $wrap-col - $eqlength ) div 2 )"/>
  <xsl:variable name="spright" select="$wrap-col - $spleft - $eqlength - $lalength"/>
  <xsl:value-of select="$newline"/>
  <xsl:value-of select="substring( '                                                                                ', 1, $spleft )"/>
  <xsl:copy-of select="$equation"/>
  <xsl:if test="$label">
    <xsl:value-of select="substring( '                                                                                ', 1, $spright )"/>
    <xsl:copy-of select="$label"/>
  </xsl:if>
  <xsl:value-of select="$newline"/>
</xsl:template>

<xsl:template name="LaTeXChar">
  <xsl:call-template name="esc">
    <xsl:with-param name="c" select='"#"'/>
    <xsl:with-param name="s">
      <xsl:call-template name="esc">
	<xsl:with-param name="c" select='"$"'/>
	<xsl:with-param name="s">
	  <xsl:call-template name="esc">
	    <xsl:with-param name="c" select='"%"'/>
	    <xsl:with-param name="s">
	      <xsl:call-template name="esc">
		<xsl:with-param name="c" select='"&amp;"'/>
		<xsl:with-param name="s">
		  <xsl:call-template name="esc">
		    <xsl:with-param name="c" select='"~"'/>
		    <xsl:with-param name="s">
		      <xsl:call-template name="esc">
			<xsl:with-param name="c" select='"_"'/>
			<xsl:with-param name="s">
			  <xsl:call-template name="esc">
			    <xsl:with-param name="c" select='"^"'/>
			    <xsl:with-param name="s">
			      <xsl:call-template name="esc">
				<xsl:with-param name="c" select='"{"'/>
				<xsl:with-param name="s">
				  <xsl:call-template name="esc">
				    <xsl:with-param name="c" select='"}"'/>
				    <xsl:with-param name="s">
				      <xsl:call-template name="esc">
					<xsl:with-param name="c" select='"\"'/>
					<xsl:with-param name="s" select='.'/>
				      </xsl:call-template>
				    </xsl:with-param>
				  </xsl:call-template>
				</xsl:with-param>
			      </xsl:call-template>
			    </xsl:with-param>
			  </xsl:call-template>
			</xsl:with-param>
		      </xsl:call-template>
		    </xsl:with-param>
		  </xsl:call-template>
		</xsl:with-param>
	      </xsl:call-template>
	    </xsl:with-param>
	  </xsl:call-template>
	</xsl:with-param>
      </xsl:call-template>
    </xsl:with-param>
  </xsl:call-template>
</xsl:template>

<xsl:template name="esc">
  <xsl:param name="s"/>
  <xsl:param name="c"/>
  <xsl:choose>
    <xsl:when test='contains($s, $c)'>
      <xsl:copy-of select='substring-before($s, $c)'/>
      <xsl:text>\</xsl:text>
      <xsl:choose>
	<xsl:when test='$c = "\"'>
	  <xsl:text>textbackslash </xsl:text>
	</xsl:when>
	<xsl:otherwise>
	  <xsl:copy-of select='$c'/>
	</xsl:otherwise>
      </xsl:choose>
      <xsl:call-template name="esc">
	<xsl:with-param name='c' select='$c'/>
	<xsl:with-param name='s' select='substring-after($s, $c)'/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:copy-of select='$s'/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template name="eqlabel">
  <xsl:param name="label"/>
  <xsl:copy-of select="count(//equation[@id=$label]/preceding::equation[@id])+1."/>
</xsl:template>

</xsl:stylesheet>
