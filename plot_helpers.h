#ifndef PLOT_HELPERS_H_
#define PLOT_HELPERS_H_

int find_last_child(int p, const std::vector<std::vector<int>>& children,
        const std::vector<int>& pdg_id) {
    int id = pdg_id[p];
    for (auto child : children[p]) {
        if (id == pdg_id[child]) {
            return find_last_child(child, children, pdg_id);
        }
    }

    return p;
}

bool is_first_parent(int p, const std::vector<std::vector<int>>& parents,
        const std::vector<int>& pdg_id) {
    auto id = pdg_id[p];
    for (auto parent : parents[p]) {
        if (id == pdg_id[parent]) {
            return false;
        }
    }
    return true;
}

TLorentzVector to_v4(int p, const std::vector<double>& pt,
        const std::vector<double>& eta, const std::vector<double>& phi,
        const std::vector<double>& m) {
    TLorentzVector v;
    v.SetPtEtaPhiM(pt[p], eta[p], phi[p], m[p]);
    return v;
}

#endif /* PLOT_HELPERS_H_ */
