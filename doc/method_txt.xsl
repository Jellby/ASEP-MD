<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:include href="txt.xsl"/>

<xsl:variable name="title">
  <xsl:copy-of select="document/title[@lang=$lang]"/>
</xsl:variable>

<xsl:template match="document">
  <xsl:value-of select="$title"/>
  <xsl:value-of select="$newline"/>
  <xsl:value-of select="$newline"/>

  <xsl:apply-templates select="paragraph[@lang=$lang]"/>
</xsl:template>

<xsl:template match="paragraph">
  <xsl:call-template name="write-block">
    <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
  </xsl:call-template>
  <xsl:value-of select="$newline"/>
  <xsl:value-of select="$newline"/>
</xsl:template>

</xsl:stylesheet>
