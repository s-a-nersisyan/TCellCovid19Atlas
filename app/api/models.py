from app import db


class HLAAlleles(db.Model):
    hla_allele = db.Column(db.Text, primary_key=True)
    gene = db.Column(db.Text, nullable=False)
    HLA_class = db.Column(db.Enum('I', 'II', name='HLA_class'), nullable=False)


class Peptides(db.Model):
    peptide = db.Column(db.Text, primary_key=True)


class PeptidesPositions(db.Model):
    peptide = db.Column(db.Text, db.ForeignKey('peptides.peptide'), primary_key=True, index=True)
    protein = db.Column(db.Text, primary_key=True, index=True)
    start = db.Column(db.SmallInteger, primary_key=True, index=True)
    end = db.Column(db.SmallInteger, primary_key=True, index=True)
    gisaid_id = db.Column(db.Text, primary_key=True, index=True)


class HLAAllelesPeptides(db.Model):
    hla_allele = db.Column(db.Text, db.ForeignKey('hla_alleles.hla_allele'), primary_key=True, index=True)
    peptide = db.Column(db.Text, db.ForeignKey('peptides.peptide'), primary_key=True, index=True)
    affinity = db.Column(db.SmallInteger, nullable=False, index=True)
