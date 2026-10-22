"""Batch encoding of gene sets with the GSFM gene-set foundation model (Clarke et al., Patterns 2026)."""
import numpy as np
import torch
from gsfm import GSFM, Vocab

PAD_IDX = 1  # GSFM special tokens: 0 = <unk>, 1 = <pad>


class Embedder:
    def __init__(self, model_id: str):
        self.vocab = Vocab.from_pretrained(model_id)
        self.model = GSFM.from_pretrained(model_id)
        self.model.eval()
        torch.set_grad_enabled(False)

    def known(self, genes):
        lk = self.vocab.lookup
        return [g for g in genes if lk.get(g, 0) > PAD_IDX]

    def coverage(self, genes) -> float:
        genes = list(genes)
        return len(self.known(genes)) / max(len(genes), 1)

    def encode(self, gene_sets, batch: int = 256) -> np.ndarray:
        out = []
        for i in range(0, len(gene_sets), batch):
            ids = [self.vocab(self.known(gs)) or [PAD_IDX] for gs in gene_sets[i:i + batch]]
            width = max(len(x) for x in ids)
            t = torch.full((len(ids), width), PAD_IDX, dtype=torch.long)
            for j, x in enumerate(ids):
                t[j, :len(x)] = torch.tensor(x, dtype=torch.long)
            out.append(self.model.encode(t).float().numpy())
        return np.vstack(out)
