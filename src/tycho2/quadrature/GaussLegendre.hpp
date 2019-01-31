template<int N>
struct GaussLegendre;

template <>
struct GaussLegendre<2> {
  static const auto& points() {
    static constexpr const std::array<double, 2> data = {
      {-0.5773502691896257, 0.5773502691896257}};
    return data;
  }
  
  static const auto& weights() {
    static constexpr const std::array<double, 2> data = {
      {1.0000000000000000, 1.0000000000000000}};
    return data;
  }
};
