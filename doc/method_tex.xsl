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

  <xsl:apply-templates select="paragraph[@lang=$lang]"/>

  <xsl:copy-of select="$enddocument"/>
</xsl:template>

<xsl:template match="paragraph">
  <xsl:apply-templates/>
  <xsl:value-of select="$newline"/>
  <xsl:value-of select="$newline"/>
</xsl:template>

</xsl:stylesheet>
