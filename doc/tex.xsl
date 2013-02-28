<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:include href="functions.xsl"/>

<xsl:output method="text" indent="no"/>

<xsl:variable name="babel">
  <xsl:choose>
    <xsl:when test="$lang='en'"><xsl:text>\usepackage[english]{babel}</xsl:text></xsl:when>
    <xsl:when test="$lang='es'"><xsl:text>\usepackage[spanish,es-noshorthands]{babel}</xsl:text></xsl:when>
  </xsl:choose>
</xsl:variable>

<xsl:variable name="begindocument">
  <xsl:text>
\documentclass[a4paper,oneside]{memoir}
\usepackage{style}
  </xsl:text>
  <xsl:value-of select="$babel"/>
  <xsl:text>
\begin{document}
  </xsl:text>
</xsl:variable>

<xsl:variable name="enddocument">
  <xsl:text>
\end{document}
  </xsl:text>
</xsl:variable>

<xsl:variable name="newline">
  <xsl:text>
</xsl:text>
</xsl:variable>

<xsl:variable name="line">
  <xsl:text>\bigskip</xsl:text>
</xsl:variable>

<xsl:template match="math">
  <xsl:text>$</xsl:text>
  <xsl:choose>
    <xsl:when test="@latex">
      <xsl:value-of select="@latex"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="."/>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>$</xsl:text>
</xsl:template>

<xsl:template match="math" mode="equation">
  <xsl:choose>
    <xsl:when test="@latex">
      <xsl:value-of select="@latex"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="."/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="exp">
  <xsl:text>\cdot 10^{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template name="tag">
  <xsl:param name="text"/>
  <xsl:text>\emph{</xsl:text><xsl:copy-of select="$text"/><xsl:text>}</xsl:text>
</xsl:template>

<xsl:template name="val">
  <xsl:param name="text"/>
  <xsl:text>\val{</xsl:text><xsl:copy-of select="$text"/><xsl:text>}</xsl:text>
</xsl:template>
<xsl:template match="val">
  <xsl:call-template name="val">
    <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
  </xsl:call-template>
</xsl:template>

<xsl:template name="var">
  <xsl:param name="text"/>
  <xsl:text>\var{</xsl:text><xsl:copy-of select="$text"/><xsl:text>}</xsl:text>
</xsl:template>
<xsl:template match="var">
  <xsl:call-template name="var">
    <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
  </xsl:call-template>
</xsl:template>

<xsl:template match="text()">
  <xsl:variable name="S"><xsl:copy-of select="."/></xsl:variable>
  <xsl:if test="string-length($S)>0">
    <xsl:call-template name="LaTeXChar"/>
  </xsl:if>
</xsl:template>

<xsl:template match="equation">
  <xsl:text>\begin{equation}</xsl:text>
  <xsl:choose>
    <xsl:when test="@id">
      <xsl:text>\label{</xsl:text>
      <xsl:value-of select="@id"/>
      <xsl:text>}</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>\notag</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:apply-templates mode="equation"/>
  <xsl:text>\end{equation}</xsl:text>
</xsl:template>

<xsl:template match="ref">
  <xsl:text>\eqref{</xsl:text>
  <xsl:value-of select="@label"/>
  <xsl:text>}</xsl:text>
</xsl:template>

</xsl:stylesheet>
