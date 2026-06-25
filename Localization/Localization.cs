namespace Crystallography;

// 260620Cl 追加: 方式② コード側多言語化ヘルパー。
// VS Designer で Localizable=false のままハードコードされたラベル (resx の ApplyResources に乗らない文字列) を、
// 実行時に現在の UI カルチャの訳へ差し替えるための N 言語版 Loc()。
// 既存の FormBeamInteraction.cs の private Loc(en, ja) を 11 言語へ一般化したもの。
// 用途・適用方針は .project-guidance/ReciPro_多言語化方針.md §3-B(方式②) を参照。
/// <summary>コード側に直書きした UI 文字列を、現在の UI カルチャの訳へ解決する localization ヘルパー。</summary>
// 260624Cl リネーム: L10n → Localization (周辺の綴り命名 CodeLocalizer / SupportedCultures に統一)
public static class Localization
{
    /// <summary>
    /// 現在の UI カルチャ (<see cref="SupportedCultures.Current"/>) に対応する訳を返す。
    /// 名前付き引数で呼ぶこと。例:
    /// <code>label.Text = Loc(en: "Space Group", ja: "空間群", de: "Raumgruppe", ...);</code>
    /// 省略・空文字の言語は英語 (<paramref name="en"/>) にフォールバックする。
    /// 引数はカルチャ名で解決するので <see cref="SupportedCultures.All"/> の並び順には依存しない
    /// (新言語を増やすときはここに 1 ケース足し、各呼び出しに該当引数を任意で足す)。
    /// </summary>
    public static string Loc(
        string en, string ja = "", string de = "", string fr = "", string es = "",
        string pt = "", string it = "", string ru = "", string zhHans = "",
        string zhHant = "", string ko = "")
    {
        var pick = SupportedCultures.Current.Name switch
        {
            "ja" => ja,
            "de" => de,
            "fr" => fr,
            "es" => es,
            "pt" => pt,
            "it" => it,
            "ru" => ru,
            "zh-Hans" => zhHans,
            "zh-Hant" => zhHant,
            "ko" => ko,
            _ => en,
        };
        return string.IsNullOrEmpty(pick) ? en : pick;
    }

    // 260621Cl 追加 (§2.5 横展開): Localizable=false のフォーム/UC が Designer.cs に英語直書きしている
    // 可視ラベル (.Text / DataGridView 列の .HeaderText / ToolStripItem の .Text) を、型名単位の中央テーブルで
    // 多言語化する。FormBase / UserControlBase の OnLoad から CodeLocalizer.Apply(this) が実行時に差し替える。
    // データ本体は自動生成の LocalizationData.cs (tools/i18n_s25_gen_data.py)。詳細は方針 §3-B/§12.7。
    /// <summary>1 コントロール (またはメニュー項目・列) の 1 プロパティぶんの 11 言語訳。</summary>
    public sealed class Entry
    {
        /// <summary>コントロール名 (Designer の Name)。フォーム自身のタイトルは "this"。</summary>
        public readonly string Ctrl;
        /// <summary>"Text" または "HeaderText"。</summary>
        public readonly string Prop;
        public readonly string En, Ja, De, Fr, Es, Pt, It, Ru, ZhHans, ZhHant, Ko;
        public Entry(string ctrl, string prop, string en, string ja, string de, string fr, string es,
                     string pt, string it, string ru, string zhHans, string zhHant, string ko)
        {
            Ctrl = ctrl; Prop = prop; En = en; Ja = ja; De = de; Fr = fr; Es = es;
            Pt = pt; It = it; Ru = ru; ZhHans = zhHans; ZhHant = zhHant; Ko = ko;
        }
        /// <summary>現在の UI カルチャに対応する訳 (<see cref="Loc"/> と同一の解決)。</summary>
        public string Resolve() => Loc(En, Ja, De, Fr, Es, Pt, It, Ru, ZhHans, ZhHant, Ko);
    }

    private static System.Collections.Generic.Dictionary<string, Entry[]> _registry;

    // 260625Cl 追加: アプリ固有の Localizable=false 訳テーブルを登録する拡張点 (multi-app 対応)。
    //   共有 LocalizationData.Populate (ReciPro/共有UC・自動生成) の後に各アプリのプロバイダを merge する。
    //   PDIndexer 等は起動時に AddProvider(PDIndexerLocalizationData.Populate) を呼び、自リポに app 固有
    //   フォームの訳を持てる (共有 Crystallography repo に app 固有フォームを混ぜない)。
    //   後方互換: プロバイダ未登録なら従来どおり LocalizationData だけ。キーは FullName 推奨 (同名衝突回避)。
    private static readonly System.Collections.Generic.List<System.Action<System.Collections.Generic.Dictionary<string, Entry[]>>> _providers = new();
    private static readonly object _sync = new();

    /// <summary>アプリ固有のローカライズ訳テーブルを登録する。フォーム生成前に1回呼ぶこと。次回 Get で再構築される。</summary>
    public static void AddProvider(System.Action<System.Collections.Generic.Dictionary<string, Entry[]>> provider)
    {
        if (provider == null) return;
        lock (_sync)
        {
            _providers.Add(provider);
            _registry = null; // 次回 Get で再構築
        }
    }

    private static void EnsureRegistry()
    {
        if (_registry != null) return;
        lock (_sync)
        {
            if (_registry != null) return;
            var reg = new System.Collections.Generic.Dictionary<string, Entry[]>(System.StringComparer.Ordinal);
            LocalizationData.Populate(reg);     // 共有 (自動生成)
            foreach (var p in _providers)
                p(reg);                          // アプリ固有 (後勝ち)
            _registry = reg;
        }
    }

    /// <summary>型名 (root.GetType().FullName または Name) に対応するローカライズ項目を返す。未登録なら null。</summary>
    public static Entry[] Get(string typeName)
    {
        EnsureRegistry();
        return typeName != null && _registry.TryGetValue(typeName, out var e) ? e : null;
    }
}
