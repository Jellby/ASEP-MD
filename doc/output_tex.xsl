<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:include href="tex.xsl"/>

<xsl:variable name="title">
  <xsl:copy-of select="document/title[@lang=$lang]"/>
</xsl:variable>

<xsl:template match="document">
  <xsl:copy-of select="$begindocument"/>

  <xsl:text>\chapter*{</xsl:text>
  <xsl:copy-of select="$title"/>
  <xsl:text>}</xsl:text>

  <xsl:apply-templates select="paragraph[@lang=$lang] | sample[@lang=$lang]"/>

  <xsl:apply-templates select="quantity"/>

  <xsl:copy-of select="$enddocument"/>
</xsl:template>

<xsl:template match="description | paragraph">
  <xsl:apply-templates/>
  <xsl:value-of select="$newline"/>
  <xsl:value-of select="$newline"/>
</xsl:template>

<xsl:template match="sample">
  <xsl:value-of select="$newline"/>
  <xsl:text>
\begin{center}
\begin{tabular}{c}
\begin{lstlisting}[style=sample]
</xsl:text>
  <xsl:copy-of select="."/>
  <xsl:value-of select="$newline"/>
  <xsl:text>
\end{lstlisting}
\end{tabular}
\end{center}
</xsl:text>
  <xsl:value-of select="$newline"/>
</xsl:template>

<xsl:template match="quantity">
  <xsl:value-of select="$line"/>
  <xsl:value-of select="$newline"/>
  <xsl:text>\oldparindent\parindent\noindent\begin{minipage}{\columnwidth}\parindent\oldparindent</xsl:text>
  <xsl:value-of select="$newline"/>

  <xsl:apply-templates select="name[@lang=$lang]"/>

  <xsl:apply-templates select="description[@lang=$lang] | description[not(@lang)]"/>

  <xsl:text>\end{minipage}</xsl:text>
  <xsl:value-of select="$newline"/>
  <xsl:value-of select="$newline"/>
</xsl:template>

<xsl:template match="name">
  <xsl:text>\noindent </xsl:text>
  <xsl:call-template name="var">
    <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
  </xsl:call-template>
  <xsl:value-of select="$newline"/>
</xsl:template>

</xsl:stylesheet>
